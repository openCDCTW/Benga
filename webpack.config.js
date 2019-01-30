module.exports = {
  module: {
    rules: [
      {
        test: /\.jsx?$/,
        exclude: /node_modules/,
        use:{
          loader: "babel-loader"
        }
      },
      {
        test:/\.(png|jpg|svg|emf|zip|pdf|tsv|newick|css|fa)$/,
        use: {
          loader: "url-loader",
        },
      }
    ]
  }
};
