#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define BUFSIZE 1024 //ファイルから読み込む一行の最大文字数
#define MAX_SEQ_NUM 30 //一つの転写因子に対して与えられる結合部位配列の最大数
#define MAX_GENE_NUM 8 //与えられるプロモータ領域の最大遺伝子数
#define THRESHOLD 2.5 //閾値

char g_motif[MAX_SEQ_NUM][BUFSIZE]; //転写因子の結合部位配列を保存する配列
int g_cols = 0; //転写因子の列数を保存する変数
double odds_score[4][MAX_SEQ_NUM]; //オッズスコア配列

struct promoter{
    char name[BUFSIZE];
    char seq[BUFSIZE];
}g_pro[MAX_GENE_NUM]; //遺伝子のプロモータ領域を保存する構造体

//グローバル変数は大文字と区別するため、全部大文字にするか、g_を先頭につけるのが一般的

//転写因子の結合部位配列を読み込む関数
int read_multi_seq(char* filename){
    int seq_num = 0;
    char buffer[BUFSIZE];
    FILE *fp = fopen(filename,"r");

    if(fp == NULL){
        printf("motif region file open error.\n");
        exit(1); //ファイルが開けなかった場合プログラムを終了
    }

    while(fscanf(fp, "%s", buffer) != EOF){ //プログラムから一行ずつ読み込む
        if(buffer[strlen(buffer)-1]=='\n'){
            buffer[strlen(buffer)-1]='\0'; //改行を切り落とす
        }
        strcpy(g_motif[seq_num],buffer); //結合部位配列を保存
        g_cols = strlen(buffer); //列数を更新
        seq_num++;
    }
    return seq_num;
}

//遺伝子のプロモータ領域を読み込む関数
int read_promoter(char *filename){
    int gene_num = 0;
    char buffer[BUFSIZE];
    FILE *fp = fopen(filename,"r");

    if(fp == NULL){
        printf("scorefile open error.\n");
        exit(1);
    }

    while(fscanf(fp, "%s", buffer) != EOF){
        if(buffer[strlen(buffer)-1]=='\n'){
            buffer[strlen(buffer)-1]='\0';
        }

        if(buffer[0]=='>'){
            strcpy(g_pro[gene_num].name,buffer+1);
        }else{
            strcpy(g_pro[gene_num].seq,buffer);
            gene_num++;
        }
    }
    return gene_num;
}

// 結合部位を探索するhit関数
double hit(char* promoter_seq, int start_pos) {
    double score = 0.0;
    
    // 各位置でのスコアを計算
    for (int i = 0; i < g_cols; i++) {
        char base = promoter_seq[start_pos + i];
        int base_index;
        
        // 塩基を数字に変換
        switch (base) {
            case 'A': base_index = 0; break;
            case 'C': base_index = 1; break;
            case 'G': base_index = 2; break;
            case 'T': base_index = 3; break;
        }
        
        // オッズスコアを加算
        score += odds_score[base_index][i];
    }
    
    return score;
}

int main(int argc, char* argv[]){
    int seq_num = read_multi_seq(argv[1]); //1番目の引数で指定した転写因子の複数の結合部位配列を読み込む

    printf("motif region:\n");
    for(int i = 0; i < seq_num; i++){
        printf("%s\n",g_motif[i]); //読み込んだ転写因子の結合部位配列を表示
    }
    printf("\n");
    printf("motif行数: %d\n", seq_num);
    printf("motif列数: %d\n", g_cols);

    int gene_num = read_promoter(argv[2]); //2番目の引数で指定した遺伝子のプロモータ領域を読み込む
/*
プロモーター領域の出力
    printf("promoter_sequence:\n");
    for(int i = 0; i < gene_num; i++){
        printf(">%s\n", g_pro[i].name); //読み込んだプロモータ領域を表示
        printf("%s\n", g_pro[i].seq);
    }
*/

//////////////////////////////////////////////////////////////////////
//頻度表の計算
//////////////////////////////////////////////////////////////////////
int freq_table[4][g_cols]; //頻度表を保存する配列
memset(freq_table, 0, sizeof(freq_table)); //初期化

// 出現頻度計算
for(int i = 0; i < seq_num; i++){
    for(int j = 0; j < g_cols; j++){
        if(g_motif[i][j] == 'A'){
            freq_table[0][j]++;
        }else if(g_motif[i][j] == 'C'){
            freq_table[1][j]++;
        }else if(g_motif[i][j] == 'G'){
            freq_table[2][j]++;
        }else if(g_motif[i][j] == 'T'){
            freq_table[3][j]++;
        }
    }
}
//頻度表の出力
printf("\n頻度表:\n");
for(int i = 0; i < 4; i++){
    for(int j = 0; j < g_cols; j++){
        printf("%4d", freq_table[i][j]);
    }
    printf("\n");
}
printf("\n");
for(int i = 0; i < g_cols; i++){
    printf("%4d", i+1);
}

// バックグラウンド出現確率
double total_bases = 7519429 + 4637676 + 4637676 + 7519429;
double background_prob_A = 7519429 / total_bases;
double background_prob_C = 4637676 / total_bases;
double background_prob_G = 4637676 / total_bases;
double background_prob_T = 7519429 / total_bases;

/*
// バックグラウンド出現確率の出力
printf("\nバックグラウンド出現確率:\n");
printf("A: %.4f\n", background_prob_A);
printf("C: %.4f\n", background_prob_C);
printf("G: %.4f\n", background_prob_G);
printf("T: %.4f\n", background_prob_T);
*/

//////////////////////////////////////////////////////////////////////
//オッズスコアの計算
//////////////////////////////////////////////////////////////////////

    for(int j = 0; j < g_cols; j++){
        double total_freq = freq_table[0][j] + freq_table[1][j] + freq_table[2][j] + freq_table[3][j] + 4;
        odds_score[0][j] = log10(((freq_table[0][j]+1) / total_freq) / background_prob_A);
        odds_score[1][j] = log10(((freq_table[1][j]+1) / total_freq) / background_prob_C);
        odds_score[2][j] = log10(((freq_table[2][j]+1) / total_freq) / background_prob_G);
        odds_score[3][j] = log10(((freq_table[3][j]+1) / total_freq) / background_prob_T);
    }
    // オッズスコアの出力
    printf("\nオッズスコア:\n");
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < g_cols; j++){
            printf("%6.2f", odds_score[i][j]);
        }
        printf("\n");
    }

    for (int gene = 0; gene < gene_num; gene++) {
        typedef struct {
            double score;
            int pos;
            char seq[BUFSIZE];
        } ScoreData;
        
        ScoreData maxScore = {-INFINITY, -1, ""};
        ScoreData aboveThreshold[BUFSIZE]; // スコアが3以上のものを保存する配列
        int count = 0; // スコアが3以上のものの数

        int seq_len = strlen(g_pro[gene].seq);
        for (int i = 0; i <= seq_len - g_cols; i++) {
            double score = hit(g_pro[gene].seq, i);
            
            // 最大スコアを更新
            if (score > maxScore.score) {
                maxScore.score = score;
                maxScore.pos = i;
                strncpy(maxScore.seq, g_pro[gene].seq + i, g_cols);
                maxScore.seq[g_cols] = '\0';
            }
            
            // スコアが閾値以上なら保存
            if (score >= THRESHOLD) {
                aboveThreshold[count].score = score;
                aboveThreshold[count].pos = i;
                strncpy(aboveThreshold[count].seq, g_pro[gene].seq + i, g_cols);
                aboveThreshold[count].seq[g_cols] = '\0';
                count++;
            }
        }

        printf("pro:%s\n", g_pro[gene].name);
        // スコアが閾値以上のものがある場合
        if (count > 0) { 
            printf("閾値%.1f以上の結合部位:\n", THRESHOLD);
            for (int i = 0; i < count; i++) {
                printf("pos:%d hit(%s)= %.2f\n", 
                    aboveThreshold[i].pos, 
                    aboveThreshold[i].seq, 
                    aboveThreshold[i].score);
            }
        // スコアが閾値以上のものがない場合
        } else { 
            printf("閾値以上がないため、最大スコア:\n");
            printf("pos:%d hit(%s)= %.2f\n", 
                maxScore.pos, 
                maxScore.seq, 
                maxScore.score);
        }
        printf("\n");
    }
    return 0;
}


