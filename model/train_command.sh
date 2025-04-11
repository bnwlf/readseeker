KMER=6
MODEL_PATH=./model/6-new-12w-0
DATA_PATH=./traindata
OUTPUT_PATH=./readseeker_finetuned/



## Try some new suggestions from  https://github.com/jerryji1993/DNABERT/issues/78
#Batchsize increased, max seq length reduced
python run_finetune.py --model_type dna --tokenizer_name=dna$KMER --model_name_or_path $MODEL_PATH --task_name dnaprom --do_train --do_eval --data_dir $DATA_PATH \
    --max_seq_length 298 --per_gpu_eval_batch_size=50 --per_gpu_train_batch_size=50 --learning_rate 2e-5 --num_train_epochs 5.0 --output_dir $OUTPUT_PATH \
    --evaluate_during_training --logging_steps 4000 --save_steps 4000 --warmup_percent 0.1 --hidden_dropout_prob 0.1 --overwrite_output --weight_decay 0.01 --n_process 11 --overwrite_cache 2>&1 |tee $OUTPUT_PATH/dnabert.log
