function c = bc_singles(param,wage,h,l,eq)

c = wage*(1-h-l)-param.pe*eq;

end