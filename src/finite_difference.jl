# Finite difference schemes
# see Driscoll 5.4 and Heath 8.6
function centered_first_deriv_wrt_pos(U, i, j, h) 
    (U[i+1,j] - U[i-1,j])/(2*h)
end 

function centered_first_deriv_wrt_time(U, i, j, h)
    (U[i,j+1] - U[i,j-1])/(2*h) 
end 

function centered_second_deriv_wrt_pos(U, i, j, h)
    (U[i+1,j] + U[i-1,j] - 2*U[i,j])/h^2
end 

function centered_second_deriv_wrt_time(U, i, j, h) 
    (U[i,j+1] + U[i,j-1] - 2*U[i,j])/h^2
end 
