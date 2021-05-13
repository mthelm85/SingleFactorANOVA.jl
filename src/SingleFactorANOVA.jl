module SingleFactorANOVA

using Distributions

export anova, AnovaResult

"""
Contains the result of calling `anova`.

# Fields
- SSB: Between treatment variation
- SSE: Error variation
- DFB: Between treatment degrees of freedom
- DFE: Error degrees of freedom
- MSB: Between treatment mean squares
- MSE: Error mean squares
- F: F-Statistic
- p: p-value for F-statistic

"""
struct AnovaResult
    SSB::Real
	SSE::Real
	DFB::Real
	DFE::Real
	MSB::Real
	MSE::Real
	F::Real
	p::Real
end

function ssb(X::AbstractArray)
    X̄ = mean(vcat(X...))
	nⱼ = [length(X[i]) for i in 1:length(X)]
	X̄ⱼ = [mean(X[i]) for i in 1:length(X)]
	sum(nⱼ .* (X̄ⱼ .- X̄).^2)
end

function sse(X::AbstractArray)
	X̄ⱼ = [mean(X[i]) for i in 1:length(X)]
	sse = []
	for j in 1:length(X)
		for i in 1:length(X[j])
			push!(sse, (X[j][i] - X̄ⱼ[j])^2)
		end
	end
	sum(sse)
end

"""
    anova(X::AbstractArray)

Perform ANOVA on `X`. `X` is an `AbstractArray` of `AbstractArray`s (typically a `Vector` of `Vectors`). Each vector
within `X` is treated as a different group and `anova` tests whether the means of the values contained in each vector 
are equal or if at least two of them are different. In other words:

    H₀: μ₁ = μ₂... = μₖ
    H₁: means are not equal

where `k` is the length of `X`.

# Example
```julia-repl
julia> result = anova([[1, 2, 5, 9], [2, 6, 4, 2, 3, 8], [15, 6, 26]])
AnovaResult(303.44230769230774, 268.25, 2, 10, 151.72115384615387, 26.825, 5.655961000788588, 0.022745050729447377)
```

To get the p-value for the F statistic:

```julia-repl
julia> result.p
0.022745050729447377
```

To view the fields of `AnovaResult`:

```julia-repl
help?> AnovaResult
search: AnovaResult

  Contains the result of calling anova.

  Fields
  ≡≡≡≡≡≡≡≡

    •  SSB: Between treatment variation
    •  SSE: Error variation
    •  DFB: Between treatment degrees of freedom
    •  DFE: Error degrees of freedom
    •  MSB: Between treatment mean squares
    •  MSE: Error mean squares
    •  F: F-Statistic
    •  p: p-value for F-statistic
```
"""
function anova(X::AbstractArray)
	N = length(vcat(X...))
	k = length(X)
	SSB = ssb(X)
	SSE = sse(X)
	DFB = k - 1
	DFE = N - k
	MSB = SSB / DFB
	MSE = SSE / DFE
	F = MSB / MSE
	fdist = FDist(DFB,DFE)
	AnovaResult(
		SSB,
		SSE,
		DFB,
		DFE,
		MSB,
		MSE,
		F,
		ccdf(fdist, F)
	)	
end

end