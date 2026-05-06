function r = wMax(m)
    ne = m.nElements;
    mm = zeros(1, ne);
    for i = 1:ne
        e = m.element(i);
        l = e.l;
        [~, v] = fminbnd(@(x) -abs(e.w(x)), 0, l);
        mm(i) = -v;
    end
    r = max(mm);
end