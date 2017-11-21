import tqdm

class TqdmBatch():
    def __init__(self, reader, batch_size, max_reads):
        self.reader = reader
        self.batch_size = batch_size
        self.progress = tqdm.tqdm(total=max_reads)

    def __iter__(self):
        for value in self.reader:
            yield value
            self.progress.update(self.batch_size)