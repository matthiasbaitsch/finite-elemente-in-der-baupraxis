out = open("_generated/basis-functions.qmd", "w")

for i in 1:length(adofs)
println(out,
"## 
Basisfunktion \$\\quad \\varphi_{$i}\$
```{julia}
plotw(
    m, 
    ei(NN, adofs[$i]), 
    w=1200, h=550,
    limits=(nothing,nothing,(-1,1))
)
\`\`\`
"
)
end


for i in 1:3
println(out,
"## 
Zufällige Kombination $i/3
```{julia}
wrand = zeros(NN)
wrand[adofs] = rand(NNa)
plotw(
    m, 
    wrand, 
    w=1200, h=550,
    limits=(nothing,nothing,(-1,1))
)
\`\`\`
"
)
end

close(out)