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

  
function spike(disc::Function, domain, p::Number)
  #δ = domain[2] - domain[1]
  y = map(disc, domain) + 0.*im
  #dy = diff(y) / δ
  dy = diff(y)
  spre = zeros(dy)
  spim = zeros(dy)
  spsq = sqrt(1 - y .^ 2)
  for i = 1:(length(dy)-1)
    spre[i+1] = spre[i] + abs(real(dy[i] / spsq[i]))
    spim[i+1] = spim[i] + imag(dy[i] / spsq[i])
  end
  return (1/p)*spre, spim
end

