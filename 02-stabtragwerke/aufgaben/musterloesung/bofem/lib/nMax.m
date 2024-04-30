function sm = nMax(m)
    ne = m.nElements;
    sms = zeros(ne, 1);
    for i = 1:ne
        sms(i) = m.element(i).N(0);
    end
    sm = max(abs(sms));
end