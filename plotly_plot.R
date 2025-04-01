data <- read.csv("data/data.csv") %>% as.matrix()
Lambda0 <- read.csv("data/lambda0.csv") %>% as.matrix()
Lambda_rotated <- read.csv("data/lambda_estimate.csv") %>% as.matrix()
truth_normal <- read.csv("data/truth_normal.csv") %>% as.matrix()

pca_vectors <- t(Lambda0)
rotated_vectors <- t(Lambda_rotated)
true_vectors <- t(truth_normal)


base_plot <- make_base_plot(data, -pca_vectors, vector_color = "black", vector_width = 6)
plotly_plot <- base_plot %>%
  add_vectors(-rotated_vectors, color = "orange", width = 6) %>%
  add_vectors(true_vectors, color = "maroon", width = 6)


  base_plot %>%
    add_vectors(-rotated_vectors, color = "orange", width = 6)
