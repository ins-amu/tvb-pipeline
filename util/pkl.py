import pickle


def read_pkl(pkl_fname):
    with open(pkl_fname, 'rb') as fd:
        obj = pickle.load(fd)
    return obj


def write_pkl(pkl_fname, obj):
    with open(pkl_fname, 'wb') as fd:
        pickle.dump(obj, fd)
