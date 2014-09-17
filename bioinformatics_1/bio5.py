import bio4

def lcs(text1, text2):
    '''
    Longest Common Substring

    !! completely wrong implementation !!

    '''
    s1, s2 = (text1, text2) if len(text1) <= len(text2) else (text2, text1)
    for n in range(len(s1), 0, -1):
        for subs in bio4.str_composition(s1, n):
            if subs in s2:
                return subs
    return ''
