function s = initRandom()
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);
end