function bits = conti2bit(pos)
    global c2bTh
    bits = pos > c2bTh;
end