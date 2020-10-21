function c = bc_married(param,wages,hf,hm,lf,lm,eq)

c = wages.f*(1-hf-lf)+wages.m*(1-hm-lm)-param.pe*eq;

end
