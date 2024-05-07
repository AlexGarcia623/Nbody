if [ -z "$1" ]; then
    # Default filename
    filename="default_simulation_params.txt"
else
    filename=$1
fi

run_dir="${filename%.txt}"
out_dir="./output/${run_dir}"

print_string="Removing all files from ${out_dir}"

length="${#print_string}"

echo "\n\n"

for (( i=0; i<length; i++ )); do
    printf "#"
done

echo "\n${print_string}"

for (( i=0; i<length; i++ )); do
    printf "#"
done

echo "\n\n"

rm "${out_dir}/particle_"*
