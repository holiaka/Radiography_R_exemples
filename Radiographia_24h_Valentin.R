# Attach the 'ForestTools' and 'raster' libraries
library("writexl")
library(ForestTools)
library(raster)

# Load a canopy height model
myRadiographia <- raster("C:\\Users\\D Holiaka\\Documents\\GitHub\\Radiography_R_exemples\\Particles_24h.tif")

# Create raster histogram before filtration
hist(myRadiographia,
     main = "Distribution of raster values",
     xlab = "Values)", ylab = "Frequency",
     col = "springgreen")


# Type main working parameters ()
working_name <- "Particles_24h"
dpi <- 508
resize_pixels <- 10
fon_level <- 2500
detect_level_for_part <- 3500
focal_filter <- 9
max_local_filter <- 20
sample_weigth_g <- 4.5

r1 <- aggregate(myRadiographia, fact = resize_pixels)

r2 <- focal(r1, w=matrix(1/focal_filter^2, nrow=focal_filter,ncol=focal_filter)) 

# Create raster histogram after filtration
hist(r2,
     main = "Distribution of raster values",
     xlab = "Values)", ylab = "Frequency",
     col = "springgreen")

# Preview
plot(r2)

# Write a crown map raster file (IF NEED)
writeRaster(r2, "C:\\Users\\D Holiaka\\Documents\\GitHub\\Radiography_R_exemples\\pre_Particles_24h.tif", dataType = "INT2U")

# Main parameters of raster layers
mean(r2)

# Main statistic of raster layers
summary(r2)

# Remove plot margins (optional)
par(mar = rep(0.5, 4))

# Plot CHM (extra optional arguments remove labels and tick marks from the plot)
plot(r2, xlab = "", ylab = "", xaxt='n', yaxt = 'n')

# Finding of local maximums
rf <- function (x){max_local_filter}
 
# Detection  of tree tops
ttops <- vwf(CHM = r2, winFun = rf, minHeight = detect_level_for_part)

# Add dominant treetops to the plot
plot(ttops, col = "blue", pch = 20, cex = 0.5, add = TRUE)

# Get the mean treetop height
min(ttops$height)
mean(ttops$height)
max(ttops$height)
# Create crown map
# IMPOTANT minHeight
crowns <- mcws(treetops = ttops, CHM = r2, minHeight = fon_level, verbose = FALSE)

# Create polygons map
# IMPOTANT minHeight
polygon <- mcws(treetops = ttops, CHM = r2, format = "polygons", minHeight = fon_level, verbose = FALSE)

# Add polygons outlines to the plot
plot(polygon, border = "red", lwd = 1.5, add = TRUE)

# Compute average area
polygon[["Area_cm2"]] <- polygon[["crownArea"]] / dpi^2 * 6.4516

# Extract raster parameters
ex_main <- extract(r2, polygon, fun=mean)
ex_median <- extract(r2, polygon, fun=median)
ex_sd <- extract(r2, polygon, fun=sd)
ex_max <- extract(r2, polygon, fun=max)
ex_sum <- extract(r2, polygon, fun=sum)

polygon[["Main"]] <- as.numeric(as.character(ex_main))
polygon[["Median"]] <- as.numeric(as.character(ex_median))
polygon[["SD"]] <- as.numeric(as.character(ex_sd))
polygon[["max"]] <- as.numeric(as.character(ex_max))
polygon[["sum"]] <- as.numeric(as.character(ex_sum))
polygon[["cor_sum"]] <- polygon[["sum"]]-(polygon[["crownArea"]]/resize_pixels^2 * fon_level)
polygon[["Activity_Bq"]] <- 3.585E-12*polygon[["cor_sum"]]^2 + 0.0000641*polygon[["cor_sum"]]
   
# Describe statistic
sp_summarise(polygon, variables = c("Area_cm2", "Main", "Median", "SD", "max", "sum", "cor_sum", "Activity_Bq"))

Pedicted_activity_Bq <- polygon$Activity_Bq
h <- hist(Pedicted_activity_Bq, xlab="Sr-90 activity Bq")
text(h$mids,h$counts,labels=h$counts, adj=c(0.5, -0.5))

# Total activity in Bq
cat("Total activity: ", sum(polygon$Activity_Bq), "Bq")

# Sr-90 activity concentration of sample
cat("Sr-90 activity concentration of sample ", sum(polygon$Activity_Bq)/sample_weigth_g*1000, "Bq/kg")

library(rgdal)

# Save a set of dominant tree tops
writeOGR(ttops, "C:\\Users\\D Holiaka\\Documents\\GitHub\\Radiography_R_exemples\\detect_paticals", working_name, driver = "ESRI Shapefile")

# Save a set of tree crown polygons
writeOGR(polygon, "C:\\Users\\D Holiaka\\Documents\\GitHub\\Radiography_R_exemples\\output_data", working_name, driver = "ESRI Shapefile")

# Add a second data set in a new worksheet
write_xlsx(polygon@data,"C:\\Users\\D Holiaka\\Documents\\GitHub\\Radiography_R_exemples\\Particles_24h.xlsx")


      