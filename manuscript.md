# A sparsity-dependent upper bound for tied-weight power autoencoders

**Michael Swissa**  
Research draft, 2 October 2026. Not submitted or published.

## Abstract

For a bias-free tied-weight autoencoder with an odd power activation and independent sparse symmetric inputs, we prove that the population reconstruction improvement is at most a constant times the hidden dimension times the sparsity probability to the reciprocal activation degree. The constant depends only on the degree and input moments. The proof combines a projection leverage inequality, two even-moment lower bounds, and a decomposition that isolates each input coordinate. Combined with the existing lower bound of Chowdhury and Weiner, this determines the optimal order throughout the sparsity range in their stated dimensional regime. No log-concavity assumption is needed for the new upper bound. Originality remains provisional pending broader expert literature review.

## Model and theorem

Fix an odd integer $m\geq3$. Let $x_i=s_i u_i$, where the $s_i$ are independent Bernoulli($p$) variables, the $u_i$ are independent with a common symmetric law $\mu$, and the two families are independent. Assume $\mu_2=\mathbb E u_i^2>0$ and $\mathbb E|u_i|^{2m}<\infty$. For even $k\leq2m$, put $\mu_k=\mathbb E u_i^k$. The distribution and $m$ are fixed when asymptotic notation is used.

For $W\in\mathbb R^{d\times n}$, set $A=W^T W$ and reconstruct $x$ as $(Ax)^{\odot m}$, where the power acts coordinatewise. Define

$$C=\mathbb E\langle x,(Ax)^{\odot m}\rangle,\qquad D=\mathbb E\|(Ax)^{\odot m}\|_2^2,\qquad I=2C-D.$$

Thus the population loss is $pn\mu_2-I$: $I$ denotes improvement over the zero reconstruction, not the loss itself. Equivalently, one may optimize over all real positive semidefinite $A$ of rank at most $d$.

**Prior ingredient.** [Scherlis et al., *Polysemanticity and Capacity in Neural Networks*](https://arxiv.org/html/2210.01892v4#A6), Appendices F and G, already establish the capacity constraint and projection leverage used below. The proposed contribution is the resulting sparsity-dependent upper bound for odd-power reconstruction, not a new rank or capacity inequality.

**Theorem.** For every $n,d$, $p\in(0,1]$, and such $A$, with $r=\operatorname{rank}A$ and $K_m=m2^{m-2}$,

$$I\leq 2K_m^2r\left(\frac{\mu_{m+1}^2}{\mu_{2m}^{(m-1)/m}\mu_2}p^{1/m}+\mu_2p\right).\tag{1}$$

In particular, $I=O(dp^{1/m})$, uniformly in $n,d,p$, since $p\leq p^{1/m}$. Neither a density nor log-concavity is assumed. At $p=0$ or $r=0$ the improvement is zero.

**Corollary (using the original paper's assumptions).** If $\mu$ is symmetric and nondegenerate with all moments finite, then for $d$ sufficiently large and $n>d^m$,

$$\sup_W I(W)=\Theta\!\left(\min\{pd^m,dp^{1/m}\}\right).\tag{2}$$

The constants and dimension threshold depend only on $m,\mu$. This uses Theorems 3.1 and 4.2 of [Chowdhury and Weiner](https://arxiv.org/html/2606.18538v1). It addresses their Section 5(1) order question for this model. We retain their stated all-moments hypothesis for the corollary and do not claim a new lower-bound construction.

The improvement matters uniformly as sparsity changes. In the strongly log-concave setting of the previous $O(d)$ bound, for example a fixed Gaussian $\mu$, at $p=d^{1-m}$ the previous upper order is $d$, whereas (1) gives order $d^{1/m}$, reducing it by a factor of order $d^{1-1/m}$. At fixed positive $p$, both upper bounds have order $d$. These statements use the unnormalized inputs $x_i=s_i u_i$ with fixed $\mu$. If inputs are instead divided by $\sqrt{p}$, the optimal improvement is divided by $p$: the corresponding matrix is rescaled by $p^{(m-1)/(2m)}$. The stated orders should not be transferred unchanged to that normalized convention.

## Proof

All expectations below are finite by the $2m$-moment assumption and finite $n$. Suppose $r>0$ and $p>0$. Let $P$ be the orthogonal projection onto the range of $A$ and write

$$\ell_i=P_{ii},\quad a_i=A_{ii},\quad q_i=\sum_j A_{ij}^2,\quad S=\sum_i a_i^m,\quad T=\sum_{i,j}A_{ij}^{2m},\quad Q=\sum_i q_i^m.$$

Positive semidefiniteness gives $a_i\geq0$. We use the previously established capacity/projection argument of Scherlis et al., Appendices F and G: because $PA=A$, Cauchy-Schwarz gives $a_i^2\leq\ell_i q_i$, while $\sum_i\ell_i=r$. In their notation the feature capacity is $a_i^2/q_i$ for $q_i>0$ (zero otherwise), so its sum is at most $r$. We reproduce the short step for completeness, without claiming novelty. Therefore

$$S^2\leq r\sum_i a_i^{2m-2}q_i\leq rT^{(m-1)/m}Q^{1/m}.\tag{3}$$

The first step applies Cauchy-Schwarz to $a_i^m\leq\sqrt{\ell_i}\,a_i^{m-1}\sqrt{q_i}$. The second applies Hölder with conjugate exponents $m/(m-1),m$ and uses $\sum_i a_i^{2m}\leq T$.

Put $D_i=\mathbb E|(Ax)_i|^{2m}$. In the multinomial expansion of this even power, independence and symmetry remove every term containing an odd coordinate exponent. Every remaining term is nonnegative, including the pure-coordinate terms. Also $\mathbb E|(Ax)_i|^2=p\mu_2 q_i$, so Jensen yields

$$D\geq p\mu_{2m}T,\qquad D_i\geq(p\mu_2q_i)^m,\qquad D\geq p^m\mu_2^mQ.\tag{4}$$

Interpolating the two global bounds with weights $(m-1)/m,1/m$ gives

$$D\geq p^{(2m-1)/m}\mu_{2m}^{(m-1)/m}\mu_2 T^{(m-1)/m}Q^{1/m}.\tag{5}$$

Combining (3) and (5), with

$$\alpha=\frac{\mu_{m+1}}{\mu_{2m}^{(m-1)/(2m)}\sqrt{\mu_2}},$$

gives

$$p\mu_{m+1}S\leq\alpha\sqrt{rD}\,p^{1/(2m)}.\tag{6}$$

Now write $(Ax)_i=a_i x_i+Z_i$, where $Z_i=\sum_{j\ne i}A_{ij}x_j$ is independent of $x_i$. Since $\mathbb E x_i=0$, $\mathbb E[x_i Z_i^m]=0$. The mean-value theorem and convexity imply, for all real $a,b$,

$$|(a+b)^m-b^m|\leq m|a|(|a|+|b|)^{m-1}\leq K_m(|a|^m+|a||b|^{m-1}).\tag{7}$$

Multiplying by $|x_i|$, taking expectations and summing yields

$$|C|\leq K_m\left(p\mu_{m+1}S+p\mu_2\sum_i a_i\mathbb E|Z_i|^{m-1}\right).\tag{8}$$

The remainder has a different bound. Conditional on $Z_i$, independence and centering give $\mathbb E[a_ix_i+Z_i\mid Z_i]=Z_i$. Conditional Jensen for the convex function $|t|^{2m}$ therefore gives $\mathbb E|Z_i|^{2m}\leq D_i$. Moment monotonicity and (4), together with $a_i^2\leq\ell_iq_i$, give

$$\mathbb E|Z_i|^{m-1}\leq D_i^{(m-1)/(2m)},\qquad a_i\leq\frac{\sqrt{\ell_i}D_i^{1/(2m)}}{\sqrt{p\mu_2}}.$$

Consequently,

$$p\mu_2\sum_i a_i\mathbb E|Z_i|^{m-1}\leq\sqrt{p\mu_2}\sum_i\sqrt{\ell_iD_i}\leq\sqrt{p\mu_2rD}.\tag{9}$$

Equations (6), (8) and (9) imply

$$|C|\leq K_m\sqrt{rD}\left(\alpha p^{1/(2m)}+\sqrt{p\mu_2}\right).\tag{10}$$

If $D>0$, completing the square gives $2C-D\leq C^2/D$, since $C^2/D-(2C-D)=(C-D)^2/D$. Squaring (10) and using $(u+v)^2\leq2(u^2+v^2)$ proves (1), without any sign assumption on $C$. If $D=0$, (4) forces every $q_i=0$, hence $A=0$ and $I=0$. This completes the proof.

For the corollary, combine (1) with the cited $O(pd^m)$ upper bound and matching $\Omega(\min\{pd^m,dp^{1/m}\})$ lower bound. The linear case is separate: for $m=1$, $I=p\mu_2\operatorname{tr}(2A-A^2)$ and $\sup I=p\mu_2\min(d,n)$, attained by a projection of that rank.

## Validation and limits

The accompanying standard-library Python program enumerates all states of small finite-support input laws and compares exact rational quantities. It tests the leverage, geometric, moment, interpolation, centering, conditional Jensen, remainder and final inequalities; verifies homogeneity by fresh enumeration; and checks positive improvements near optimal scalar rescaling. Input powers are $3,5,7$, with $n\leq4$, deterministic and seeded random Gram matrices, and three symmetric laws, including an atom at zero. The machine-readable results give all counts and timings. These finite checks are falsification attempts, not a proof or evidence of training behavior.

A fresh AI mathematical reviewer independently rederived the theorem and found no gap; a further review pass audited the written notation and tests. These are not external human peer review or formal proof verification. The exact run passed 576 cases, with 74,986 repeated moment-enumeration state visits and 311,525 counted checks, including 405 strictly positive scalar-rescaled improvements. State visits are not unique inputs and exclude separate linear enumerations. The final-bound and combined-correlation checks all take the nonpositive-delta branch because the theorem's constant is loose at these small dimensions; they add little stress. The intermediate inequalities provide the substantive falsification checks. The three controls are direct counterexamples to false strengthenings or violated centering, rather than mutations run through the whole harness. The result concerns asymptotic order, not sharp finite-dimensional constants.

The proof requires independent coordinates and independent Bernoulli masks, a real positive semidefinite tied-weight matrix, no bias, and an odd integer power activation. It does not establish the broader conditional-expectation result in the paper's Section 5(2), allow arbitrary dependent masks, prove training convergence, or imply performance guarantees for language models, general autoencoders, or deep networks. No model was trained.

## Sources, priority and disclosure

The primary reference is Mriganka Basu Roy Chowdhury and Eric McLaughlin Weiner, *Effects of sparsity and superposition on loss in simple autoencoders*, [arXiv:2606.18538v1](https://arxiv.org/html/2606.18538v1), 16 June 2026, Sections 2-5 and Lemma A.4. The [arXiv version history](https://arxiv.org/abs/2606.18538) showed only v1 when checked on 2 October 2026; the [author's research page](https://mriganka.xyz/research) linked that preprint. Their lower bound and model supply the context for (2); the elementary upper-bound proof above is the proposed contribution.

Targeted public searches for the paper title, identifier, authors, and sparsity-dependent power-autoencoder bounds did not identify an earlier proof of (1). A recent citing paper, [*Interference Beyond Geometry in Concept Extraction*](https://arxiv.org/html/2609.35351v1), was checked for relevance; its concept-extraction framework did not provide this bound. These checks do not certify worldwide priority. A novelty claim should remain provisional until a specialist conducts a broader citation and literature review.

An important direct antecedent is Adam Scherlis, Kshitij Sachan, Adam S. Jermyn, Joe Benton and Buck Shlegeris, [*Polysemanticity and Capacity in Neural Networks*, arXiv:2210.01892v4](https://arxiv.org/html/2210.01892v4#A6), 25 March 2025 revision (first posted in 2022). Their capacity/rank constraint is an existing ingredient of this proof. Their scalar quadratic-regression model differs from the odd-power vector reconstruction considered here; we did not identify bound (1) in that work.

Classical frame-energy inequalities are not claimed as new. Relevant background includes Martin Ehler and Kasso A. Okoudjou, [*Minimization of the Probabilistic p-frame Potential*](https://arxiv.org/abs/1101.0140), and Shayne Waldron, [*A variational characterisation of projective spherical designs over the quaternions*](https://www.math.auckland.ac.nz/~waldron/Preprints/Quaternion-designs/quaternion-designs.pdf), Theorem 4.1 and Corollary 5.1. No frame-design appendix or unverified construction is needed for the main theorem.

**AI assistance.** This draft, proof audit, literature searches and test code were developed with OpenAI Codex AI assistance, including separate AI reviewer sessions. Michael Swissa is the supplied author identity; no affiliation or credentials are asserted, and no claim is made that he personally derived or checked the mathematics. This is an AI-assisted research draft, without external peer review. Independent human expert review is recommended before relying on the result. Only public sources and synthetic inputs were used for the research; no private-project content, employer data, paid API, author contact, public repository or publication was used.
