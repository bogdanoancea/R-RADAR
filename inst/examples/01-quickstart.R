library(RADAR)

dataset <- generate(4, length = 128, seed = 7)
components <- RADAR::decompose(dataset$X, model = "standard")
print(components)
print(evaluate_components(dataset, components))

output <- file.path(tempdir(), "radar-components.rds")
RADAR::write(components, output)
restored <- load_components(output)
stopifnot(isTRUE(all.equal(restored$peaks, components$peaks)))
plot(components, index = 1)
