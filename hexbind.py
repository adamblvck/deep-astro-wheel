hexagram_bin = np.floor(neutron_stream) # rounded up downwards

# map bin number onto 'hexagram' (neutron stream is sequential order, hexagram is King Wen Sequence)
strong = np.array(iching_map)
previous_shape = hexagram_bin.to_numpy().shape
mapped = strong[flat]
hexagram = mapped.reshape(previous_shape)



