# Assortative-mating

Some R code to simulate data assuming assortative mating (phenotypic correlation between partners) arisen through 3 non-mutually exclusive processes.

1. individual-specific trait values (i) of males and female are correlated between partners

2. there is a common environment that affects both partners and hence creates correlated common environmental effects

3. there are (reciprocal) indirect effects between the partners. 

There are several hierarchical levels affecting the assortative mating correlation. Alomng the lines of earlier work partitioning the P matrix. 
Question is whether one can construct the assortative mating correlation on the level of individuals 
Suppose a simple model, where there is a common environment (c.e.)

z_mal = i_mal + c.e_mal + epsilon_mal

z_fem = i_fem + c.e_fem + epsilon_fem


WORKING SCENARIO 1
TRUE ASSORTATIVE MATING. The common environmental effect is not correlated for the partners 

cov(c.e_mal, c.e_fem) = 0
cov(i_mal, i_fem) > 0

Assuming var(i_mal), cov(i_mal,i_fem), var(i_fem) are from a multivariate normal distribution

There are each year n males paired to n females. They form a pair for 3 years. Then all the females die and new females are mated to. After this pairs change, they are together for 3 years again. Then all the males die and the females remade and breed for 3 years. 
Total of 9 years with 9 * 1000 = 9000 observation of which 3000 with an animal monogamous and 6000 with 2 partners for both sexes.
1.	Draw i_mal and i_fem from MVN
2.	for each of 3 year, draw c.e and epsilon 
3.	add up to get trait values
4.	Remate
5.	step 2 and 3
6.	remate
7.	step 2 and 3

The above is done by AssortativeMatingSim.R

Second, analyse the results
a.	from the raw data calcuklate the “traditional” assortative mating
b.	use the individual.specific means of the data to calculate assortative mating
c.	estimate the mixed model
c(trait.mal, trait.fem) ~ Trait + at(trait,1).malID + at(trait,2).femID  + MVN(res)
where the assorative mating is the correlation in P (or a separate model of that)
d.	use the model BLUPs to estimate assortative mating

Step 3.	construct a model in JAGS or STAN where the cov(i.mal,ifem) can be specified.

i.	draw from variance males
ii.	use cov to estimate females
iii.	for ind with second partner, do this again (both sexes)
cov approach described in D&R.
some careful data sorting is needed with help vectors
the major problem is to simulate re-mating of individuals that have survived. that is, find from a pool of individuals with ind values, the ones that correlate. still, for first tritls the above described “death of partner” approach should work

Data structure
m1	f1	(3 rep measures per pair)
m2	f2	(1000 instead of 2)
m1	f3
m2	f4
m3	f3
m4	f4
males mated once = c(m3,m4); mated multiple times = c(m1,m2)
females mated once= c(f1,f2) ; mated multiple times = c(f3,f4)
1.	draw from var males with multiple times
2.	use cov to estimate their females’ value (single + multiple)
3.	use the values for the females mated multiple times to estimate the value of the males mated once with them.

??? can any structure be reduced to this?? simple to add individuals only mated to each other (strictly monogamous pairs) as these are given by the covar directly.
? mated with more than 2 partners, does it then get impossibly messy?
 the real challenge is when individuals remate with each other (re-partnering pn the pool of existing individuals. how then to decide which trait values to draw????? We need a function that allows to do this
