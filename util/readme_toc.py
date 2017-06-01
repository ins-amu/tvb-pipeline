def print_readme_toc():

    with open('README.md', 'r') as fd:
        lines = fd.readlines()

    tocs = []
    for line in lines:
        line = line.strip()
        if line[:2] == '##':
            parts = line.split()
            level = len(parts[0])
            title = ' '.join(parts[1:])
            slug = '-'.join([p.lower() for p in parts[1:]])
            tocs.append((level, title, slug))

    print('**Table of contents**\n')
    for level, title, slug in tocs:
        space = '  ' * (level - 2)
        fmt = '%s- [%s](#%s)'
        fmt %= space, title, slug
        print(fmt)

if __name__ == '__main__':
    print_readme_toc()
    # TODO replace toc directly..
