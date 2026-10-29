library(radar.raman)

dataset <- generate(4, length = 128, seed = 7)
components <- radar.raman::decompose(dataset$X, model = "standard")
print(components)
print(evaluate_components(dataset, components))

output <- file.path(tempdir(), "radar-components.rds")
radar.raman::write(components, output)
restored <- load_components(output)
stopifnot(isTRUE(all.equal(restored$peaks, components$peaks)))
plot(components, index = 1)
