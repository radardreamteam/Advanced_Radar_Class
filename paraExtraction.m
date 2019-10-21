function [tx2tg,rx2tg] = paraExtraction(Ddifference, tx2rx, angle)

tx2tg = ((tx2rx+Ddifference)^2+ tx2rx^2- 2*tx2rx*Ddifference*cos(angle))...
/(2*(tx2rx+Ddifference-tx2rx*cos(angle)));
rx2tg =  ((tx2rx+Ddifference)^2- tx2rx^2)...
        /(2*(tx2rx+Ddifference-tx2rx*cos(angle)));
end
    