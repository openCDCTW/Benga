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
        test:/\.(png|jpg|svg)$/,
        use: {
          loader: "url-loader"
        }
      }
    ]
  }
};
