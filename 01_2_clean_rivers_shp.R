
# clean the rivers data a bit so that nodes match up correctly 
library(sf)

rivers <- read_sf(file.path(Sys.getenv("DATA_DIR"), "SWORD_v16_gpkg", "na_sword_reaches_v16.gpkg"))


rg <- st_geometry(rivers)

# first and second points overlap for all vectors -- mostly
# loop through rows to make sure
for (row in 1:nrow(rivers)){
  id <- rivers[row,]$reach_id
  upstream <- as.numeric(strsplit(rivers[row,]$rch_id_up, " ")[[1]])
  # want: last(down) == first(up)
  for (up in upstream){
    # first chec to see if it's even in the dataset
    if (up %in% rivers$reach_id) {
      whichup <- which(rivers$reach_id == up)
      if ( 
        # don't change if last(down) == first(up) already
        !all( rg[[row]][nrow(rg[[row]]), ] == rg[[whichup]][1, ] )
      ){
        
        # if (last-1)(down) == first(up) trim last on down
        if (
          all( rg[[row]][nrow(rg[[row]])-1, ] == rg[[whichup]][1, ] )
        ) {
          # trim last row
          rg[[row]] <- st_linestring(rg[[row]][1:(nrow(rg[[row]])-1),])
        } else if (
          # if (last-2)(down) == first(up) trim last two on down
          all( rg[[row]][nrow(rg[[row]])-2, ] == rg[[whichup]][1, ] )
        ) {
          # trim last two rows
          rg[[row]] <- st_linestring(rg[[row]][1:(nrow(rg[[row]])-2),])
        } else if (
          # last(down) == last(up)
          all( rg[[row]][nrow(rg[[row]]), ] == rg[[whichup]][nrow(rg[[whichup]]), ] )
        ) {
          # reverse order of geometry for upstream reach
          rg[[whichup]] <- st_linestring(rg[[whichup]][nrow(rg[[whichup]]):1,])
        } else if (
          # (last-1)(down) == last(up)
          all( rg[[row]][nrow(rg[[row]])-1, ] == rg[[whichup]][nrow(rg[[whichup]]), ] )
        ) {
          # trim last row
          rg[[row]] <- st_linestring(rg[[row]][1:(nrow(rg[[row]])-1),])
          # reverse
          rg[[whichup]] <- st_linestring(rg[[whichup]][nrow(rg[[whichup]]):1,])
        } else if (
          # (last-2)(down) == last(up)
          all( rg[[row]][nrow(rg[[row]])-2, ] == rg[[whichup]][nrow(rg[[whichup]]), ] )
        ) {
          # trim last row
          rg[[row]] <- st_linestring(rg[[row]][1:(nrow(rg[[row]])-2),])
          # reverse
          rg[[whichup]] <- st_linestring(rg[[whichup]][nrow(rg[[whichup]]):1,])
        } else {
          print("could not match up for: ")
          print(paste("down: ", id))
          print(paste("up:", up))
        }
      }
    }
    
  }
  
}

# reverse direction so it goes downstream
rg <- rg |>
  lapply(\(g) st_linestring(g[nrow(g):1,])) |>
  st_sfc(crs = st_crs(rivers))

st_geometry(rivers) <- rg

# Save to disk (in lat/lon for geojson)
rivers |>
  st_transform(crs= "EPSG:4326") |>
  st_write(here::here("data/rivers.geojson"))






