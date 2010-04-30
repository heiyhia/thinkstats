scores = [55, 66, 77, 88, 99]
your_score = 88

def PercentileRank(scores, your_score):
    count = 0
    for score in scores:
        if score <= your_score:
            count += 1

    percentile_rank = 100.0 * count / len(scores)
    return percentile_rank


for score in scores:
    print score, PercentileRank(scores, score)


def Percentile(scores, percentile_rank):
    scores.sort()
    for score in scores:
        if PercentileRank(scores, score) >= percentile_rank:
            return score

def Percentile2(scores, percentile_rank):
    scores.sort()
    index = percentile_rank * (len(scores)-1) / 100
    return scores[index]

for percentile_rank in [0, 20, 25, 40, 50, 60, 75, 80, 100]:
    print percentile_rank, Percentile2(scores, percentile_rank)
