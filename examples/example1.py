import maxwell_slice as ms

if __name__ == '__main__':
    data, field = ms.sphere_scat()
    print(data.shape)
    print(field.shape)