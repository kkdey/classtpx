
## limSolve tutorial

E <- matrix(nrow = 3, ncol = 3, data = c(3, 1, 2, 2, 0, 0, 1, 0, 2))
F <- c(2, 1, 8)
solve(E, F)




E2 <- rbind(E, E[1, ] + E[2, ])
F2 <- c(F, F[1] + F[2])
solve(E2,F2)  ## gives error
Solve(E2, F2) ## gives the correct output: same as above for the determined system

## over-determined system

## Overdetermined systems contain more independent equations than unknowns
## the least squares solution can be extracted using the function lsei()

##  for equalities case

A <- matrix(nrow = 4, ncol = 3, data = c(3, 1, 2, 0, 2, 0, 0, 1, 1, 0, 2, 0))
B <- c(2, 1, 8, 3)
lsei(A = A, B = B, fulloutput = TRUE, verbose = FALSE)

## for equalities + inequalities

G <- matrix(nrow = 2, ncol = 3, byrow = TRUE, data = c(-1, 2, 0, 1, 0, -1))
H <- c(-3, 2)
lsei(A = A, B = B, G = G, H = H, verbose = FALSE)

## Under-determined system

## Underdetermined systems contain less independent equations than unknowns
## we aim for the most parsimonious solution

## equalities case

E <- matrix(nrow = 2, ncol = 3, data = c(3, 1, 2, 0, 1, 0))
F <- c(2, 1)
Solve(E, F)
ldei(E = E, F = F)$X

## equalities + inequalities case

E <- matrix(ncol = 4, byrow = TRUE, data = c(3, 2, 1, 4, 1, 1, 1, 1))
F <- c(2, 2)
G <-matrix(ncol = 4, byrow = TRUE, data = c(2, 1, 1, 1, -1, 3, 2, 1, -1, 0, 1, 0))
H <- c(-1, 2, 1)
ldei(E, F, G = G, H = H)$X
