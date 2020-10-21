function u = obj_married(param,ctrl,pwf)

h = h_married(param,ctrl(3),ctrl(4));

g = home_prod(param,h,ctrl(7));
            
u = u_married(param,ctrl(1),ctrl(2),ctrl(5),ctrl(6),g,pwf);

end

