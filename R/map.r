library(sf)
library(mapsf)
library(spData) # For the us_states dataset
library(airportr)
library(bezier)
library(h3jsr)
library(dplyr)
data("na_states")

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

#' Merge the states data with a CSV file containing additional information (e.g., days spent in each state) based on the GEOID column.
#' @param states An `sf` object containing the state geometries and a GEOID column.
#' @param days_csv An optional path to a CSV file containing additional state data with a GEOID column. If NULL, an internal example CSV will be used.
#' @param crs An optional coordinate reference system to transform the states data to after merging. If NULL, the original CRS will be retained.
#' @return An `sf` object with the merged data.
#' @export
merge_states <- function(states, by="iso_3166_2", days_csv=NULL, crs=NULL) {
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

  # Merge the states data with the CSV data based on the specified column
  states_data = read.csv(file_path, stringsAsFactors=FALSE)
  states = merge(states, states_data, by="iso_3166_2", all.x=TRUE)
 
  # Optionally transform the CRS of the states data
  if (!is.null(crs)) {
    states <- st_transform(states, crs)
  }

  return(states)
}

#' Draw a choropleth map of states colored by a specified variable
#' @param states An `sf` object containing the state geometries and a GEOID column.
#' @param var The name of the variable in the states data to use for coloring.
#' @param pal A vector of colors to use for the choropleth.
#' @param breaks A numeric vector of break points for the choropleth.
#' @param add A logical value indicating whether to add to an existing map. Default is FALSE.
#' @export
draw_states <- function(states, var, pal, breaks, add = FALSE) {
  mf_map(x = states, type = "base", col = NA, border = NA, add = add)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#9abac6", border = NA) 
  
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

  arcs <- st_transform(st_as_sf(flights, geometry = st_sfc(curves, crs = 4326)), st_crs(states))
  arcs <- st_wrap_dateline(arcs, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
  airports = st_transform(st_as_sf(airports, coords = c("Longitude", "Latitude"), crs = 4326), st_crs(states))

  return(list(airports = airports, arcs = arcs))
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

  bbox <- st_as_sfc(st_bbox(states)) 
  bbox <- st_segmentize(bbox, units::set_units(10, km))
  bbox <- st_transform(bbox, 4326) # H3 operates in WGS84 (EPSG:4326)
  bbox <- st_wrap_dateline(bbox, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=10"))
  bbox <- st_buffer(bbox, dist = 250000)
  
  hex_ids <- polygon_to_cells(bbox, res = res, simple = TRUE)

  h3_hexes <- cell_to_polygon(hex_ids, simple = FALSE)
  h3_hexes <- st_transform(h3_hexes, st_crs(states))

  # find all the visited hexes by checking for intersections with the provided sf objects (e.g., airports, roads)
  list_of_logical_vectors <- lapply(spatial_layers, function(layer) {
      # Ensure CRS matches h3_hexes and drop Z if needed
      layer_clean <- sf::st_zm(sf::st_transform(layer, st_crs(h3_hexes)))
  
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
# Gather and transform the data
#
data("us_states")
png(
    filename = "map.png",
    width = 2560,
    height = 1920,
    res = 150,
    type = "cairo")

mf_init(x = st_transform(us_states, 5070), expandBB = c(0.15, 0, 0, 0)) 

states_csv = system.file("extdata", "states.csv", package = "travelmapr")
flights_csv = system.file("extdata", "flights.csv", package = "travelmapr")
roadtrips_dir = system.file("extdata", "roadtrips", package = "travelmapr")

states <- merge_states(na_states, days_csv = states_csv, crs=5070)
flights <- get_flights(states, flights_csv = flights_csv, k=0.45, k_noise=0.1)
roads <- get_road_trips(states, roadtrips_dir)

hexes <- get_hexes(states, spatial_layers = list(flights$airports, roads), res=3)

#
# Draw the map
#
pal <- c("white","#fefeb5", "#fee391", "#fec44f", "#fe9929", "#d95f0e", "#993404", "#662506" ) 
breaks <- c(0, 0.041, 1, 7, 30, 365, 3650, 36500)

draw_states(states, var = "days", pal = pal, breaks = breaks, add = TRUE)

mf_map(roads, col = "#0000ff88", lwd=0.25, add = TRUE )
mf_map(hexes, col="#FFFFFF78", border="#00000014", lwd=0.5, add = TRUE)
mf_map(flights$airports, col = "black", pch = 18, cex = 0.5, add = TRUE)
mf_map(flights$arcs, col = "#000000aa", lwd = 0.25, add = TRUE)

# Hawaii inset
hawaii <- states[states$iso_3166_2 == "US-HI",]
st_agr(hawaii) <- "constant"
hawaii <- st_crop(st_transform(hawaii, 4326),
                  xmin = -161, xmax = -154, ymin = 18, ymax = 23)
hawaii <- st_transform(hawaii, 3759)

mf_inset_on(fig = c(0.36, 0.55, 0.005, 0.125))
draw_states(hawaii, var = "days", pal = pal, breaks = breaks)
hi_hexes <- get_hexes(hawaii, spatial_layers = list(flights$airports, roads), res=3)
box(col = "black", lwd = 1) 

mf_map(hi_hexes, col="#FFFFFF98", border="#00000014", lwd=0.5, add = TRUE)
mf_map(st_transform(roads, st_crs(hawaii)), col = "#0000ff88", lwd=0.25, add = TRUE )
mf_map(st_transform(flights$airports, st_crs(hawaii)), col = "black", pch = 18, cex = 0.5, add = TRUE)
mf_map(st_transform(flights$arcs, st_crs(hawaii)), col = "#000000aa", lwd = 0.25, add = TRUE)
mf_inset_off()

# Alaska inset
alaska <- st_transform(states[states$iso_3166_2 == "US-AK",], 3338)
alaska_region <- st_transform(states[states$iso_3166_2 %in% c("US-AK","CA-YT","CA-BC","CA-NT"),], 3338)

mf_inset_on(x = alaska, pos = "bottomleft", cex = 0.35) 
mf_init(x = alaska)
ak_hexes <- get_hexes(alaska, spatial_layers = list(flights$airports, roads), res=3)
draw_states(alaska_region, var = "days", pal = pal, breaks = breaks, add=TRUE)
box(col = "black", lwd = 1) 
mf_map(ak_hexes, col="#FFFFFF98", border="#00000014", lwd=0.5, add = TRUE)
mf_map(st_transform(roads, st_crs(alaska)), col = "#0000ff88", lwd=0.25, add = TRUE )
mf_map(st_transform(flights$airports, st_crs(alaska)), col = "black", pch = 18, cex = 0.5, add = TRUE)
mf_map(st_transform(flights$arcs, st_crs(alaska)), col = "#000000aa", lwd = 0.25, add = TRUE)
mf_inset_off()

# Legend 
mf_inset_on(fig = c(0.54, 1.1, 0.0, 0.13)) 
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
rect(0, 0, 1, 1, col = "#cadae1", border = "black")

mf_legend(
  type = "choro", pos = c(0.025, 0.40), val = c("","hours", "days", "weeks", "months", "years", "decades", ""),
  pal = pal, title = "", horiz = TRUE, val_cex = 1.0, title_cex = 1.0, size=1.95
)
text(0.845, 0.86, "April 2026", pos = 4, cex = 1.15)
text(0.895, 0.0825, "by TravelMapR", pos = 4, cex = 0.5)


x <- seq(0.24, 0.34, length.out = 20)
y <- 0.6 + 0.1 * sin(seq(0, pi, length.out = 20))
lines(x, y, col = "black", lwd = 0.75)
text(0.345, 0.6, "By Air", pos = 4, cex = 0.9)

lines(c(0.48, 0.58), c(0.6, 0.6), col = "darkblue", lwd = 0.5, lty = 1)
text(0.585, 0.6, "By Land", pos = 4, cex = 0.9)

lines(c(0.745, 0.845), c(0.6, 0.6), col = "darkred", lwd = 0.5, lty = 1)
text(0.855, 0.6, "By Sea", pos = 4, cex = 0.9)

points(0.075, 0.6, col = "black", pch = 18, cex = 1.0)
text(0.0775, 0.6, "Airports", pos = 4, cex = 0.9)

mf_inset_off()

dev.off()
graphics.off()
