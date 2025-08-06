mkdir -p "simple_model_results/timings/"


for sample in data/*;
do
samplen=$(basename $sample | cut -f 2 -d ".");
ref=$(basename $sample | cut -f 1 -d ".");

echo $ref
echo $samplen
echo $sample

"time" -f "SimpleBenchmarkModel,$samplen,$ref,%U,%S,%e" -o "simple_model_results/timings/"$samplen"_"$ref".txt" python benchmark_simple_benchmark_model.py simple_benchmark_model_saved_final $sample $ref $samplen;
done

cat simple_model_results/timings/*.txt >  simple_model_results/timings.txt
