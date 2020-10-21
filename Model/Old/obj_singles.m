function u = obj_singles(param,wage,controls)

c = bc_singles(param,wage,controls(1),controls(2),controls(3));
g = home_prod(param,controls(1),controls(3));

u = u_singles(param,c,controls(2),g);

end


