library(sf)
library(mapsf)
library(spData) # For the us_states dataset
library(airportr)
library(bezier)
library(h3jsr)
library(dplyr)

#' Create a Bezier arc (i.e., stylistic flight path) between two points, with an optional curvature parameter k.
#'
#' @param p1 A numeric vector of length 2 representing the first point (longitude, latitude).
#' @param p2 A numeric vector of length 2 representing the second point (longitude, latitude).
#' @param n The number of points to generate along the curve.
#' @param k The curvature parameter.
#' @param k_noise The amount of noise to add to the curvature parameter for visual variety.
#' @return An `sf` linestring representing the Bezier curve.
#' @export
make_bezier_arc <- function(p1, p2, n = 50, k = 0.5, k_noise = 0.1) {

  # Add some noise to k for visual variety
  k <- runif(1, k-k_noise, k+k_noise) 

  # Handle the Date Line crossing
  # If the points are more than 180 degrees apart horizontally,
  # shift one of them temporarily so the math calculates the short path.
  lon1 <- p1[1]
  lon2 <- p2[1]
  
  if (abs(lon1 - lon2) > 180) {
    if (lon1 < 0) lon1 <- lon1 + 360
    if (lon2 < 0) lon2 <- lon2 + 360
  }

  p1_mod <- c(lon1, p1[2])
  p2_mod <- c(lon2, p2[2])

  # Create a control point for the Bezier curve
  mid <- (p1_mod + p2_mod) / 2
  v <- p2_mod - p1_mod

  if (p2_mod[1] < p1_mod[1]) {
    k <- -k # arcs up
  }

  perpendicular <- c(-v[2], v[1]) * k 
  perpendicular[2] <- pmin(perpendicular[2], 10.0) 
  control_pt <- mid + perpendicular

  # Generate the curve points
  pts <- bezier(t = seq(0, 1, length = n), 
                p = rbind(p1_mod, control_pt, p2_mod))

  # Bring any shifted longitudes back into the standard [-180, 180] range
  pts[, 1] <- ifelse(pts[, 1] > 180, pts[, 1] - 360, pts[, 1])
  pts[, 1] <- ifelse(pts[, 1] < -180, pts[, 1] + 360, pts[, 1])

  return(st_linestring(pts))
}

#' Draw a choropleth map of states colored by a specified variable, with an optional CSV file for state-level data.
#' @param states An `sf` object containing the state geometries and a GEOID column.
#' @param var The name of the variable in the states data to use for coloring.
#' @param pal A vector of colors to use for the choropleth.
#' @param breaks A numeric vector of break points for the choropleth.
#' @param days_csv An optional path to a CSV file containing state-level data with a GEOID column. If NULL, an internal example CSV will be used.
#' @return A merged `sf` object containing the state geometries joined with the data from the CSV file.
#' @export
draw_states <- function(states, days_csv = NULL, var, pal, breaks) {
  # Read the CSV file, either from the internal example or from the user-provided path
  if (is.null(days_csv)) {
    # Pull the internal example file
    file_path <- system.file("extdata", "states.csv", package = "travelmapr")
    message("Using internal package example data...")
  } else {
    # Use the path provided by the user
    file_path <- days_csv
  }

  if (!file.exists(file_path)) {
    stop(paste0("File not found: '", file_path, "'. Please check the path and try again."), 
         call. = FALSE)
  }

  # Merge the states data with the CSV data based on the GEOID column
  states$GEOID = as.numeric(states$GEOID)
  states_data = read.csv(file_path, stringsAsFactors=FALSE)
  states = merge(states, states_data, by="GEOID", all.x=TRUE)

  # Draw the map
  mf_map(x = states, type = "base", col = NA, border = NA)

  mf_shadow(
    x = states,
    col = "grey50",
    cex = 0.5,
    add = TRUE
  )

  mf_map(
    x = states,
    var = var,
    type = "choro",
    pal = pal,
    breaks = breaks,
    nbreaks = length(breaks) - 1,
    leg_pos = NA,
    add = TRUE
  )

  return(states)
}

#' Compute flight paths between airports based on a CSV file of flight data, 
#' with optional parameters curvature, and noise.
#' @param flights_csv An optional path to a CSV file containing flight data with two columns for origin and destination IATA codes. If NULL, an internal example CSV will be used.
#' @param k The curvature parameter for the Bezier arcs representing flight paths.
#' @param k_noise The amount of noise to add to the curvature parameter for visual variety in the flight paths.
#' @return A list wtih `sf` objects containing the airport locations as points and the flight paths (arcs) as linestrings.
#' @export
get_flights <- function(states, flights_csv=NULL, k=0.45, k_noise=0.1) {
  if (is.null(flights_csv)) {
    # Pull the internal example file
    file_path <- system.file("extdata", "flights.csv", package = "travelmapr")
    message("Using internal package example data...")
  } else {
    # Use the path provided by the user
    file_path <- flights_csv
  }

  if (!file.exists(file_path)) {
    stop(paste0("File not found: '", file_path, "'. Please check the path and try again."), 
         call. = FALSE)
  }

  flights = read.csv(file_path, header=FALSE, stringsAsFactors = FALSE) 

  airports = data.frame("IATA" = unique(c(flights[,1], flights[,2])))
  airports = merge(airports, airportr::airports[c('IATA', 'Latitude', 'Longitude')], by="IATA", all.x=TRUE) 

  curves <- lapply(1:nrow(flights), function(i) {
    make_bezier_arc(
      c(airports$Longitude[airports$IATA == flights[i, 1]], airports$Latitude[airports$IATA == flights[i, 1]]),
      c(airports$Longitude[airports$IATA == flights[i, 2]], airports$Latitude[airports$IATA == flights[i, 2]]),
      k = k, k_noise = k_noise
    )
  })

  arcs <- st_as_sf(flights, geometry = st_sfc(curves, crs = 4326))
  arcs <- st_wrap_dateline(arcs, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))

  return(list(airports = st_as_sf(airports, coords = c("Longitude", "Latitude"), crs = st_crs(states)), arcs = arcs))
}


#' Create road trips on the map based on KML files in a specified directory
#' @param states An `sf` object containing the state geometries and a GEOID column, used for coordinate reference.
#' @param roadtrip_dir An optional path to a directory containing KML files of road trips. If NULL, an internal example directory will be used.
#' @return An `sf` object containing the road trip paths as linestrings.
#' @export
get_road_trips <- function(states, roadtrip_dir=NULL) {
  if (is.null(roadtrip_dir)) {
    roadtrip_dir <- system.file("extdata", "roadtrips", package = "travelmapr")
  }

  if (!dir.exists(roadtrip_dir)) {
    stop(paste0("Directory not found: '", roadtrip_dir, "'. Please check the path and try again."), 
         call. = FALSE)
  }

  kml_files <- list.files(roadtrip_dir, pattern = "\\.kml$", full.names = TRUE)
  
  roads_list <- lapply(kml_files, function(f) {
    lyrs <- st_layers(f)$name
    file_data <- lapply(lyrs, function(l) {
      st_read(f, layer = l, quiet = TRUE)
    })
    do.call(rbind, file_data)
  })
  
  roads <- bind_rows(roads_list)
  roads <- roads[st_is(roads, "LINESTRING"), ]
  roads <- st_transform(roads, st_crs(states))
  
  return(roads)
}

#' Get hexagonal grid cells over the map, keeping only those that do not intersect with any of the provided spatial features (e.g., airports, roads).
#' @param states An `sf` object containing the state geometries, used to define the area for the hexagonal grid.
#' @param sf_list A list of `sf` objects representing spatial features (e.g., airports, roads) that should be checked for intersections with the hexagonal grid cells.
#' @param res The resolution parameter for the H3 hexagonal grid (higher values create smaller hexagons).
#' @return An `sf` object containing the hexagonal grid cells that do not intersect with any of the provided spatial features.
#' @export
get_hexes <- function(states, spatial_layers, res=3) {

  bbox <- st_buffer(st_transform(st_as_sfc(st_bbox(states)), 5070), 50000) # buffer in meters to ensure coverage beyond the borders of the states
  bbox <- st_transform(bbox, 4326) # H3 operates in WGS84 (EPSG:4326)
  hex_ids <- polygon_to_cells(bbox, res = res, simple = TRUE)
  h3_hexes <- cell_to_polygon(hex_ids, simple = FALSE)
  h3_hexes <- st_transform(h3_hexes, st_crs(states))
  h3_hexes <- h3_hexes[st_intersects(h3_hexes, st_union(states), sparse = FALSE), ] # only keep hexes that intersect with the states

  # find all the visited hexes by checking for intersections with the provided sf objects (e.g., airports, roads)
  list_of_logical_vectors <- lapply(spatial_layers, function(layer) {
      # Ensure CRS matches h3_hexes (4326) and drop Z if needed
      layer_clean <- sf::st_zm(sf::st_transform(layer, st_crs(states)))
  
      # Check for intersections
      inter <- sf::st_intersects(h3_hexes, layer_clean, sparse = FALSE)
  
      # Return a single TRUE/FALSE per hex
      apply(inter, 1, any)
  })

  has_intersections <- Reduce(`|`, list_of_logical_vectors)
  h3_hexes <- h3_hexes[!has_intersections, ]

  return(h3_hexes)
}

#
# Draw the map
#
data("us_states")

states_csv = system.file("extdata", "states.csv", package = "travelmapr")
flights_csv = system.file("extdata", "flights.csv", package = "travelmapr")
roadtrips_dir = system.file("extdata", "roadtrips", package = "travelmapr")

pal <- c("white","#fefeb5", "#fee391", "#fec44f", "#fe9929", "#d95f0e", "#993404", "#662506" ) 
breaks <- c(0, 0.041, 1, 7, 30, 365, 3650, 36500)

png(
    filename = "map.png",
    width = 2560,
    height = 1640,
    res = 150,
    type = "cairo")

states <- draw_states(us_states, days_csv = states_csv, var = "days", pal = pal, breaks = breaks)
flights <- get_flights(states, flights_csv = flights_csv, k=0.45, k_noise=0.1)
roads <- get_road_trips(states, roadtrips_dir)
hexes <- get_hexes(states, spatial_layers = list(flights$airports, roads), res=3)
mf_map(roads, col = "#0000ff88", lwd=0.25, add = TRUE )
mf_map(hexes, col="#FFFFFF98", border="#00000014", lwd=0.5, add = TRUE)
mf_map(flights$airports, add = TRUE, col = "black", pch = 18, cex = 0.5)
mf_map(flights$arcs, col = "black", lwd = 0.25, add = TRUE)

dev.off()
graphics.off()

