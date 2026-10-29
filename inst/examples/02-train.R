library(radar.raman)

sets <- train_test_datasets(4, 2, length = 64, seed = 3)
model <- new_model("standard", scaler = 2, dim_red_size = 16, seed = 3)
model <- train(sets[[1]], model = model, epochs = 1, batch_size = 2,
               validation_data = sets[[2]], learning_rate = 1e-5, verbose = 1)
print(model$history)
path <- file.path(tempdir(), "radar-trained.rds")
save_model(model, path)
restored <- load_model(path)
print(evaluate(restored, sets[[2]]))
