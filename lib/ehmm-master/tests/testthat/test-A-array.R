context("A: Array routines")

#R EQUIVALENTS TO THE COMPILED C++ FUNCTIONS
#returns a list with two elements:
#1: values, sorted unique values of the counts vector
#2: map, reference to the values vector such that:
# all(values[map+1] == counts)
their_mapToUnique <- function(counts){
	o <- order(counts)
	uval <- rle(counts[o])
	values <- uval$values
	uval$values <- 1:length(values)
	map <- integer(length(counts))
	map[o] <- inverse.rle(uval) - 1
	storage.mode(map) <- "integer"
	storage.mode(values) <- "integer"
	list(values=values, map=map)
}

#START TESTING

test_that("Array routines work",{
	#we don't set the seed here
	counts <- exampleData(1000, indep=T)
	
	#mapToUnique
	their_ucss <- apply(counts, 1, their_mapToUnique)
	my_ucss <- apply(counts, 1, mapToUnique)
	expect_identical(their_ucss, my_ucss)
	
	#colSumsInt
	their_cs <- colSums(counts); storage.mode(their_cs) <- "integer"
	my_cs <- colSumsInt(counts)
	expect_identical(my_cs, their_cs)
	
	#colSumsDouble and rowSumsDouble
	rcounts <- matrix(rnorm(1e4), nrow=10)
	their_cs <- colSums(rcounts)
	my_cs <- colSumsDouble(rcounts)
	expect_equal(my_cs, their_cs)
	
	their_rs <- rowSums(rcounts)
	my_rs <- rowSumsDouble(rcounts)
	expect_equal(my_rs, their_rs)
})
