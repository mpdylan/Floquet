module PlotFloq

function findbandUC(disc::Function, domain)
  counter = 1
  if abs(disc(domain[1])) < 2 
    while true
      if abs(disc(domain[counter])) < 2
        counter += 1
      else
        break
      end
    end
    while true
      if abs(disc(domain[counter])) > 2
        counter += 1
      else
        break
      end
    end
  else
    while true
      if abs(disc(domain[counter])) > 2
        counter += 1
      else
        break
      end
    end
  end
  counter
end

function spikeUC(disc::Function, domain)
  border = findbandUC(disc, domain)
  n = length(domain)
  y = zeros(Complex{Float64}, n)
  z = disc(domain[border])
  for i=border:n-1
    w = disc(domain[i+1])
    δ = w - z
    sq = sqrt(4. + 0.0im - z^2)
    y[i + 1] = y[i] + sign(δ) * real(δ / sq) - im * sign(z) * imag(δ / sq)
    z = w
  end
  if border > 1
    w = disc(domain[1])
    δ = w - z
    sq = sqrt(4. + 0.0im - z^2)
    y[1] = y[n] + sign(δ) * real(δ / sq) - im * sign(z) * imag(δ / sq)
    z = w
    for i=1:border - 1
      w = disc(domain[i+1])
      δ = w - z
      sq = sqrt(4. + 0.0im - z^2)
      y[i + 1] = y[i] + sign(δ) * real(δ / sq) - im * sign(z) * imag(δ / sq)
      z = w
    end
  end
  y
end

function adspikeUC(disc::Function, tol::Float64)
  if tol <= 0
    error("Tolerance must be positive")
  end
  currentx = 0
  dx = tol
  mintol = tol
  y = 0.
  #disc(x) = real(Floquet.discUC(α, exp(im * x)))
  function nextstep()
    dy = abs(disc(currentx + dx) - disc(currentx))
    dx *= tol/dy
    if dx < mintol
      mintol = dx
    end
    #print(dx)
    #y = disc(currentx + dx)
    #currentx += dx
  end
  if abs(disc(0)) < 2                #First we seek the first lower band edge
    while true                       #If we started inside a band we have to proceed through the next gap
      nextstep()
      y = disc(currentx + dx)
      currentx += dx
      #print(currentx, "\n")
      if abs(y) > 2
        break
      end
    end
  end
  while true                         #The next part is to locate the lower band edge. If we didn't start in a band we can begin here
      nextstep()
      y = disc(currentx + dx)
      currentx += dx
      #print(currentx, "\n")
      if abs(y) < 2
        break
      end
  end
  print(currentx)
                                     #Now we compute the spike map
  start = currentx
  stepest = round(2π / mintol) 
  output1 = zeros(Complex{Float64}, Int64(round(stepest*(1 - currentx / (2π)))))
  output0 = zeros(Complex{Float64}, Int64(round(stepest*currentx / (2π))))
  y = disc(currentx)
  y2 = disc(currentx + dx) 
  dy = y2 - y
  i = 1
  while currentx < 2π
    sq = sqrt(4. + 0.0im - y^2)
    if i < length(output1)
      output1[i+1] = output1[i] + sign(dy) * real(dy / sq) - im * sign(y) * imag(dy / sq)
      i += 1
    else
      nextout = output1[i] + sign(dy) * real(dy / sq) - im * sign(y) * imag(dy / sq)
      output1 = vcat(output1, nextout)
      i += 1
    end
    nextstep()
    y = y2
    currentx += dx
    y2 = disc(currentx)
  end
  if i <= length(output1)
    resize!(output1, i-1)
  end
  currentx = 0
  i = 1
#  sq = sqrt(4. + 0.0im - y^2)
#  output0[1] = output1[end] + sign(dy) * real(dy / sq) - im * sign(y) * imag(dy / sq)
  while currentx < start
      sq = sqrt(4. + 0.0im - y^2)
    if i < length(output0)
      output0[i+1] = output0[i] + sign(dy) * real(dy / sq) - im * sign(y) * imag(dy / sq)
      i += 1
    else
      nextout = output0[i] + sign(dy) * real(dy / sq) - im * sign(y) * imag(dy / sq)
      output0 = vcat(output0, nextout)
      i += 1
    end
    nextstep()
    y = y2
    currentx += dx
    y2 = disc(currentx)
  end
  if i <= length(output0)
    resize!(output0, i-1)
  end
  vcat(output0, output1)
end


  
#function spike(disc::Function, domain, p::Number)
  #δ = domain[2] - domain[1]
#  y = map(disc, domain) + 0.*im
  #dy = diff(y) / δ
#  dy = diff(y)
#  spre = zeros(dy)
#  spim = zeros(dy)
#  spsq = sqrt(1 - y .^ 2)
#  for i = 1:(length(dy)-1)
#    spre[i+1] = spre[i] + abs(real(dy[i] / spsq[i]))
#    spim[i+1] = spim[i] + imag(dy[i] / spsq[i])
#  end
#  return (1/p)*spre, spim
#end

end

