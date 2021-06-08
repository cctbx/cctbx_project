

def conversion():
    input = open("C:\\Users\\Ken\\Desktop\\Richardson\\suitename-duke\\janesviews.txt")

    lines = input.readlines()
    trimmed = [line[2:-4] for line in lines]
    text = "\n".join(trimmed)
    print(text)


# conversion()

def conversion2():
    input = open("helptext.txt")
    lines = input.readlines()
    trimmed = []
    for line in lines:
        begin = line.find('"') + 1
        end = line.find('"', begin)
        core = line[begin:end]
        trimmed.append(core)
    text = "\n".join(trimmed)
    print(text)

conversion2()

