build_B = function(vec, N){
  #get size of the vector
  dim = length(vec)
  #Partition the vector of presence absence it in rows of length N
  #vec_mat = matrix(vec, nrow = N, byrow = T)
  #Preallocate matrix B
  B = matrix(0, N, dim)
  for (i in seq(N)){
    for (k in seq(N-1)){
      #In row i, I can only modify the elements in columns of the form j = 1+N*i
      j = i + N*(k - 1)
      #If row k has a 1 in element i, then B[i, j] = 1, where j = 1+N*i
      vec_slice = vec[(N*(k-1)+1):(N*k)]
      if (vec_slice[i] == 1){
        B[i, j] = 1
      }
    }
  }
  #If row k has a 1 in element i, then B[i, j] = 1, where j = 1+N*i
  #Else, B[i, j] = 0
  return(B)
}

B = build_B(vec, N)
