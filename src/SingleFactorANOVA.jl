module SingleFactorANOVA

using Distributions

export anova, AnovaResult, tukey_kramer, TukeyKramerResult

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

"""
Contains the result of calling `tukey-kramer`.

# Fields
- mean_difference: Pairwise differences between means
- q_crit: Pairwise critical q values
- significant: Boolean values indicating whether pairwise mean differences are significant
"""
struct TukeyKramerResult
	mean_differences::NamedTuple
	q_crit::NamedTuple
	significant::NamedTuple
end

"""
    tukey_kramer(X::AbstractArray, anova_result::AnovaResult, α=0.05)

Perform Tukey-Kramer test on `X`. `X` is the same `AbstractArray` that you passed to the `anova` function. This test
is only performed if you have rejected the null hypothesis after having called `anova`. It returns a `TukeyKramerResult` which
shows which pairs of means are significantly different from one another.

# Example
```julia-repl
julia> X = [[1, 2, 5, 9], [2, 6, 4, 2, 3, 8], [15, 6, 26]]
3-element Vector{Vector{Int64}}:
 [1, 2, 5, 9]
 [2, 6, 4, 2, 3, 8]
 [15, 6, 26]

julia> result = anova(X)
AnovaResult(303.44230769230774, 268.25, 2, 10, 151.72115384615387, 26.825, 5.655961000788588, 0.022745050729447377)

julia> julia> tk_result = tukey_kramer(X, result, 0.05)
TukeyKramerResult((|x1 - x2| = 0.08333333333333304, |x1 - x3| = 11.416666666666666, |x2 - x3| = 11.5), (|x1 - x2| = 9.164737679909036, |x1 - x3| = 10.843863861104227, |x2 - x3| = 10.03946712180748), (|x1 - x2| = false, |x1 - x3| = true, |x2 - x3| = true))
```

To view the fields of `TukeyKramerResult`:

```julia-repl
help?> TukeyKramerResult
search: TukeyKramerResult

  Contains the result of calling tukey-kramer.

  Fields
  ≡≡≡≡≡≡≡≡

    •  mean_difference: Pairwise differences between means

    •  q_crit: Pairwise critical q values

    •  significant: Boolean values indicating whether pairwise mean differences are significant
```

To just see the `significant` field:

```julia-repl
julia> tk_result.significant
(|x1 - x2| = false, |x1 - x3| = true, |x2 - x3| = true)
```
"""
function tukey_kramer(X, anova_result, α=0.05)
	x̄ = [abs(mean(X[i]) - mean(X[j])) for i in 1:length(X) for j in i+1:length(X)]	
	d = StudentizedRange(anova_result.DFE, length(X))
	q_crit = [quantile(d, 1 - α) * sqrt((anova_result.MSE/2)*((1/length(X[i]))+(1/length(X[j])))) for i in 1:length(X) for j in i+1:length(X)]
	comparisons = Tuple(Symbol("|x$i - x$j|") for i in 1:length(X) for j in i+1:length(X))
	TukeyKramerResult(
		NamedTuple{comparisons}(x̄),
		NamedTuple{comparisons}(q_crit),
		NamedTuple{comparisons}(x̄ .> q_crit)
	)
end

end