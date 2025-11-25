function p_pos_norm = PositivePressures(p)
% Transforming for pressure to positive pressure. Since compressor only
% can provide positive pressure, need to push the pressure to the positive
% scale.
% 

minP = min(p);
maxP = max(p);

p_pos = p-minP;
p_pos_norm = p_pos/max(p_pos) * maxP;

end