# Clustering for coculture phase 
using RDatasets, Clustering, Plots

# loading the dataset
iris = dataset("datasets", "iris");

# features for clustering
features = collect(Matrix(iris[:, 1:4])');

# result after running K-means for the 3 clusters
result = kmeans(features, 3);

# plotting the result
scatter(iris.PetalLength, iris.PetalWidth,
		marker_z = result.assignments,
		color =:lightrainbow, legend = false)

# saving the result in PNG form
# savefig("D:\\iris.png")
