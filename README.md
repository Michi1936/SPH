# SPH
粒子法（SPH法のソースコード）
SPH.h:パラメータ、関数の定義
numbers.h:流体粒子、壁粒子の数を記述
SPH.c:密度、圧力などの各物理量を計算
setting.c:流体粒子、壁粒子を配置
Bucket.c:計算を高速にするため、計算領域を格子状に分割し近傍の粒子との未計算を行うようにする
main.c:上述の関数を呼び出す
particle.py画像データを流体、壁粒子の初期位置の座標に変換
anime.txt:出力されたデータをgnuplot上で再生するためのコード
fluid.txt, wall.txt:流体粒子、壁粒子の初期位置

コンパイル:
gcc main.c SPH.c setting.c Bucket.c -lm
