cmd.bg_color("white")
                                               
load ./tst_00_answer.pdb, answer
cmd.hide("everything","poor")
cmd.show("sticks"    ,"poor")
cmd.show("spheres"   ,"poor")
color red, answer

load ./tst_00_refined_all_states.pdb, rsr
cmd.hide("everything","rsr")
cmd.show("sticks"    ,"rsr")
cmd.show("spheres"   ,"rsr")
color black, rsr


load ./map.ccp4, map1, 1 , ccp4
isomesh mesh1, map1, 1.5, (all), 0, 1, 3.0
cmd.color("blue","mesh1")
                                                 # 
                                                 # 
                                                 # 
set stick_radius,0.04                              
util.cbaw("pmod")                                  
set_color cyan,[ 0.38, 0.78, 0.80]                 
set_color blue,[ 0.11, 0.25, 0.88]                 
set_color red,[ 1.00, 0.13, 0.13]                  
#color cyan, element h & mod6                      
#color grey, element c & mod6                      
#color red, all # for the final picture only       
#color black, name ca | name n | name c | name O   
set sphere_scale,0.075                             
                                                   
set mesh_width,0.35                                
set mesh_radius,0.001 #0.012                       
                                                   
set direct,1.
set orthoscopic,on
set ray_trace_fog_start,0.5
util.performance(0)
util.ray_shadows('light')
cmd.space('cmyk')
set antialias,on
                                                  #
#set_view (\
#    -0.755569279,    0.136784866,    0.640635848,\
#    -0.646187246,   -0.316175640,   -0.694611311,\
#     0.107539862,   -0.938792348,    0.327278197,\
#     0.000000000,    0.000000000,  -54.740055084,\
#     8.479999542,    9.727499962,    9.920499802,\
#    42.835460663,   66.644645691,   20.000000000 )
#rebuild
#
