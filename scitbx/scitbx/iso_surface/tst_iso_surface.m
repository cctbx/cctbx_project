f[x_,y_,z_]:=Sin[x*y+y*z+z*x]
f[x,y,z]
D[f[x,y,z],{x,2}]//Simplify
D[f[x,y,z],{y,2}]//Simplify
D[f[x,y,z],{z,2}]//Simplify
D[f[x,y,z],x,y]//Simplify
D[f[x,y,z],y,z]//Simplify
D[f[x,y,z],z,x]//Simplify