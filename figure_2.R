# Plots for constants of motion in two dimensions
# written by Abigail D'Ovidio Long

# Parameters for fast transition rate
a_12 = 0.45
lambda = 0.4

# fixed constants 
z_1 = 0.4
z_2 = 0.6
z_3 = 0.8
z_4 = 0.2


#initiate values
y_01   <- 0:1
mat <- cbind(A = 1/lambda*(z_1*(a_12+lambda) - a_12*y_01), C = 1/lambda*(z_2*(a_12+lambda) - a_12*y_01), D = 1/lambda*(z_3*(a_12+lambda) - a_12*y_01), C = 1/lambda*(z_4*(a_12+lambda) - a_12*y_01), B = y_01)

# plot 
cols <- c("red", "red", "red", "red", 'black')
par(pty = "s")
matplot(y_01, mat, type = "l", lty = 1, lwd = 2, col = cols, 
        xlim = c(0,1), ylim = c(0,1),  xaxs = "i",
        yaxs = "i", main = "Constants of Motion, n = 2",
        xlab = expression(y["0,1"]), ylab = expression(y["0,2"]))
arrows(x0 = 0.35, y0 = 0.03125, x1 = .3, y1 = .0875, # values calculated from lines
       col = "red",          
       lwd = 2,               
       length = 0.1,         
       code = 2)              # sets head at start
arrows(x0 = 0.05, y0 = .36875, x1 = .1, y1 = .3125, 
       col = "red",          
       lwd = 2,              
       length = 0.1,         
       code = 2)              
arrows(x0 = 0.55, y0 = .23125, x1 = .5, y1 = .2875, 
       col = "red",          
       lwd = 2,               
       length = 0.1,         
       code = 2)             
arrows(x0 = 0.25, y0 = .56875, x1 = .3, y1 = .5125, 
       col = "red",         
       lwd = 2,              
       length = 0.1,        
       code = 2)              
arrows(x0 = 0.75, y0 = .43125, x1 = .7, y1 = .4875, 
       col = "red",          
       lwd = 2,               
       length = 0.1,        
       code = 2)              
arrows(x0 = .45, y0 = .76875, x1 = .5, y1 = .7125, 
       col = "red",          
       lwd = 2,               
       length = 0.1,         
       code = 2)              
arrows(x0 = .95, y0 = .63125, x1 = .9, y1 = .6875, 
       col = "red",         
       lwd = 2,               
       length = 0.1,         
       code = 2)              
arrows(x0 = 0.65, y0 = .96875, x1 = .7, y1 = .9125, 
       col = "red",          
       lwd = 2,               
       length = 0.1,         
       code = 2)              
