function sm = sigmaMax(m)
    ne = m.nElements;
    sms = zeros(ne, 1);
    for i = 1:ne
        sms(i) = m.element(i).sigmaMax();
    end
    sm = max(sms);
end
