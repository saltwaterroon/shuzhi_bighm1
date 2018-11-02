//
//  main.cpp
//  MyFirstOpenCV
//
//  Created by h-xiao16 on 2018/10/24.
//  Copyright © 2018年 h-xiao16. All rights reserved.
//


#include <opencv2/highgui.hpp>
#include <opencv2/opencv.hpp>
#include"mypoint.h"
using namespace cv;
using namespace std;
const int n=68;

vector<mypoint>src_pt;//原图片的点
vector<mypoint>tg_pt;//目标图像的点
vector<mypoint>src_ft;//原图片的特征点
vector<mypoint>tg_ft;//目标图像特征点
//定义基函数u(r)
static double tps_base_func(double r)
{
    if ( r == 0.0 )
        return 0.0;
    else
        return r*r * log(r*r);
}
double *ALU(double a[n+3][n+3], double b[n+3])
{
    double l[n+3][n+3] = { 0 };
    double u[n+3][n+3] = { 0 };
    int i, r, k;
    //进行U的第一行的赋值
    for (i = 0; i<n+3; i++)
    {
        u[0][i] = a[0][i];
    }
    
    //进行L的第一列的赋值
    for (i = 1; i<n+3; i++)
    {
        l[i][0] = a[i][0] / u[0][0];
    }
    
    //计算U的剩下的行数和L的剩下的列数
    for (r = 1; r<n+3; r++)
    {
        for (i = r; i <n+3; i++)
        {
            double sum1 = 0;
            for (k = 0; k < r; k++)
            {
                sum1 += l[r][k] * u[k][i];
                //cout << "" << r << "" << sum1 << endl;
            }
            u[r][i] = a[r][i] - sum1;
        }
        
        
        if(r!=n+3)
            for(i=r+1;i<n+3;i++)
            {
                double sum2 = 0;
                for (k = 0; k<r; k++)
                {
                    sum2 += l[i][k] * u[k][r];
                }
                l[i][r] = (a[i][r] - sum2) / u[r][r];
            }
        
    }
    
    double y[n+3] = { 0 };
    y[0] = b[0];
    for (i = 1; i<n+3; i++)
    {
        double sum3 = 0;
        for (k = 0; k<i; k++)
            sum3 += l[i][k] * y[k];
        y[i] = b[i] - sum3;
    }
    
    double * x= new double [n+3];
    x[n+3]={0};
    x[n+2] = y[n+2] / u[n+2][n+2];
    for (i = n+1; i >= 0; i--)
    {
        double sum4 = 0;
        for (k = i + 1; k<n+3; k++)
            sum4 += u[i][k] * x[k];
        x[i] = (y[i] - sum4) / u[i][i];
    }
    for (i = 0; i<n+3; i++)
        cout << "x[" << i + 1 << "]=" << x[i] << endl;
    return x;
}

double * gaussin(double **a, double *b)
{
//   double** a=(double**)malloc((n+3)*sizeof(double*));
//  double*  b=(double*)malloc((n+3)*sizeof(double));
//  //  double *mtx_v2_change=(double*)malloc((n+3)*sizeof(double));
//    // w1=ALU(mtx_l,mtx_v1);
//    //w2=ALU(mtx_l,mtx_v2);
//    for(int i=0;i<n+3;i++)
//    {
//        a[i]=(double*)malloc((n+3)*sizeof(double));
//       
//    }
//    a=a1;
//    b=b1;
    for(int i=0;i<n+3;i++)
    {
        for(int j=0;j<n+3;j++)
            cout<<a[i][j]<<" ";
        cout<<endl;
        
    }
    //判断能否用高斯消元法，如果矩阵主对角线上有0元素存在是不能用的
//    for (int i = 0; i < n+3; i++)
//        if (a[i][i] == 0)
//        {
//                        cout<<"no."<<i<<endl;
//            cout << "can't use gaussin meathod" << endl;
//            return 0;
//
//        }
//    
    int i, j, k;
    double c[n+3];    //存储初等行变换的系数，用于行的相减
    //消元的整个过程如下，总共n-1次消元过程。
    for (k = 0; k < n +3- 1; k++)
    {
        //求出第K次初等行变换的系数
        for (i = k + 1; i < n+3; i++)
            c[i] = a[i][k] / a[k][k];
        
        //第K次的消元计算
        for (i = k + 1; i < n+3; i++)
        {
            for (j = 0; j < n+3; j++)
            {
                a[i][j] = a[i][j] - c[i] * a[k][j];
            }
            b[i] = b[i] - c[i] * b[k];
        }
    }
    
    //解的存储数组
    double *x=new double[n+3];
    //先计算出最后一个未知数；
    x[n+3 - 1] = b[n+3 - 1] / a[n+3 - 1][n +3- 1];
    //求出每个未知数的值
    for (i = n+3 - 2; i >= 0; i--)
    {
        double sum = 0;
        for (j = i + 1; j < n+3; j++)
        {
            sum += a[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / a[i][i];
    }
    
    cout << " the solution of the equations is:" << endl;
    cout << endl;
    for (i = 0; i < n+3; i++)
        cout <<"x"<<i+1<<"="<< x[i] << endl;

    return x;
    
}
double* Gause(double** a, int n)
{
    int i, j, k;
    int rank, columm;
    double temp,l,s;
    double*x = new double[n];
    for (i = 0; i <= n - 2; i++)
    {
        rank=i;
        columm=i;
        for(j=i+1;j<=n-1;j++)                     //选主元
            if (a[j][i] > a[i][i])
            {
                rank = j;
                columm = i;
            }
        for (k = 0; k <= n; k++)                //主元行变换
        {
            temp = a[i][k];
            a[i][k] = a[rank][k];
            a[rank][k] = temp;
        }
        for (j = i + 1; j <= n - 1; j++)         //解线性方程
        {
            l = a[j][i] / a[i][i];
            for (k = i; k <= n; k++)
                a[j][k] = a[j][ k] - l * a[i][k];
        }
    }
    x[n - 1] = a[n - 1][n] / a[n - 1][ n - 1];       //回代方程求解x
    for (i = n - 2; i >= 0; i--)
    {
        s = 0;
        for (j = i + 1; j <= n - 1; j++)
            s = s + a[i][j] * x[j];
        x[i] = (a[i][ n] - s) / a[i][i];
        
    }
    cout<<"解为";
    for(int i=0;i<n;i++)
    {
        cout<<i<<":"<<x[i]<<endl;
    }
    return x;
}

double ** tps_result()
{
    
   // unsigned n = tg_ft.size();
    
    double mtx_l[n+3][n+3];
   // mtx_l[71][71]={0};
//    double mtx_v[n+3][2];
//    mtx_v[n+3][2]={0};
    double mtx_k[n][n];
 //   mtx_k[n][n]={0};
    double mtx_v1[n+3];
 //   mtx_v1[n+3]={0};
    double mtx_v2[n+3];
 //   mtx_v2[n+3]={0};
    //K
    for ( unsigned i=0; i<n; ++i )
    {
        for ( unsigned j=i+1; j<n; ++j )
        {
            mypoint pt_i = tg_ft[i];
            mypoint pt_j = tg_ft[j];
            //pt_i.y = pt_j.y = 0;
            double uij = (pt_i - pt_j).distance();
            mtx_l[i][j] = mtx_l[j][i] =mtx_k[i][j]= mtx_k[j][i] =tps_base_func(uij);
        }
    }
    
    for ( unsigned i=0; i<n; ++i )
    {
        //对角线
        mtx_l[i][i] = mtx_k[i][i] =0;
        //p
        mtx_l[i][n] = 1.0;
        mtx_l[i][n+1]= tg_ft[i].x;
        mtx_l[i][n+2]= tg_ft[i].y;
        mtx_l[n][i] = 1.0;
        mtx_l[n+1][i]= tg_ft[i].x;
        mtx_l[n+2][i]= tg_ft[i].y;
    }
    
    
    // O (3 x 3, lower right)
    for ( unsigned i=n; i<n+3; ++i )
        for ( int j=n; j<n+3; ++j )
            mtx_l[i][j] = 0.0;

   
    // 方程右边 V
    for ( unsigned i=0; i<n; ++i )
    {
        mtx_v1[i]= src_ft[i].x;
        mtx_v2[i]= src_ft[i].y;
    }
    
    for ( unsigned i=n; i<n+3; ++i )
    {
        mtx_v1[i]= 0;
        mtx_v2[i]= 0;
    }
    double **mtx_l_change=(double**)malloc((n+3)*sizeof(double*));
    double *mtx_v1_change=(double*)malloc((n+3)*sizeof(double));
    double *mtx_v2_change=(double*)malloc((n+3)*sizeof(double));
    for(int i=0;i<3;i++)
    {
        mtx_l_change[i]=(double*)malloc((n+3)*sizeof(double));
        mtx_l_change[i]=mtx_l[n+i];
        mtx_v1_change[i]=mtx_v1[n+i];
        mtx_v2_change[i]=mtx_v2[n+i];
    }
    for(int i=3;i<n+3;i++)
    {
        mtx_l_change[i]=(double*)malloc((n+3)*sizeof(double));
        mtx_l_change[i]=mtx_l[i-3];
        mtx_v1_change[i]=mtx_v1[i-3];
        mtx_v2_change[i]=mtx_v2[i-3];
    }
    double **mtx_increase=(double**)malloc((n+3)*sizeof(double*));
    
    for(int i=0;i<n+3;i++)
    {
        mtx_increase[i]=(double*)malloc((n+4)*sizeof(double));
        for(int j=0;j<n+3;j++)
            mtx_increase[i][j]=mtx_l_change[i][j];
        mtx_increase[i][n+3]=mtx_v1_change[i];
    }
   
    double **mtx_increase2=(double**)malloc((n+3)*sizeof(double*));
    
    for(int i=0;i<n+3;i++)
    {
        mtx_increase2[i]=(double*)malloc((n+4)*sizeof(double));
        for(int j=0;j<n+3;j++)
            mtx_increase2[i][j]=mtx_l_change[i][j];
        mtx_increase2[i][n+3]=mtx_v2_change[i];
    }

//    for(int i=0;i<n+2;i++)
//    {
//        mtx_l_change[i]=(double*)malloc((n+3)*sizeof(double));
//        mtx_l_change[i]=mtx_l[i+1];
//        mtx_v1_change[i]=mtx_v1[i+1];
//        mtx_v2_change[i]=mtx_v2[i+1];
//    }
//    mtx_l_change[n+2]=(double*)malloc((n+3)*sizeof(double));
//    mtx_l_change[n+2]=mtx_l[0];
//    mtx_v1_change[n+2]=mtx_v1[0];
//    mtx_v2_change[n+2]=mtx_v2[0];
    ofstream outfile3;
    outfile3.open("/Users/h-xiao16/Desktop/数值分析与算法/mtxl_change.txt"); //mtxl.txt是存放数据的文件名
    
    for ( unsigned i=0; i<n+3; ++i )
    {
        for ( unsigned j=0; j<n+4; ++j )
        {
            printf("%f",mtx_increase[i][j]);
            printf(" ");
            if(outfile3.is_open())
            {
                outfile3<<mtx_increase[i][j]<<" ";    //message是程序中处理的数据
                
                
            }
            
            
            
        }
        printf("\n");
        outfile3<<endl;
    }
    
    outfile3.close();
//    //mtx_l换行
//    for(int i=0;i<n+3;i++)
//    {
//        int a;
//        a=mtx_l[0][i];
//        mtx_l[0][i]=mtx_l[1][i];
//        mtx_l[1][i]=a;
//    }

    int b;
    b=mtx_v1[0];
    mtx_v1[0]= mtx_v1[1];
    mtx_v1[1]= b;
    b=mtx_v2[0];
    mtx_v2[0]= mtx_v2[1];
    mtx_v2[1]= b;
    //输出k
    ofstream outfile;
    outfile.open("/Users/h-xiao16/Desktop/数值分析与算法/mtxl.txt"); //mtxl.txt是存放数据的文件名
  
    for ( unsigned i=0; i<n+3; ++i )
    {
        for ( unsigned j=0; j<n+3; ++j )
        {
            printf("%f",mtx_l[i][j]);
            printf(" ");
            if(outfile.is_open())
            {
                outfile<<mtx_l[i][j]<<" ";    //message是程序中处理的数据
              
                
            }
            
        
            
        }
        printf("\n");
        outfile<<endl;
    }
    
    outfile.close();
//    for ( unsigned i=0; i<n; ++i )
//    {
//        for ( unsigned j=0; j<n; ++j )
//        {
//            printf("%f",mtx_k[i][j]);
//            printf(" ");
//            
//        }
//        printf("\n");
//    }
    ofstream outfilev1;
    outfilev1.open("/Users/h-xiao16/Desktop/数值分析与算法/mtxv1.txt"); //mtxl.txt是存放数据的文件名

    printf("v is");
      for ( unsigned i=0; i<n+3; ++i )
      {
          printf("%f",mtx_v1[i]);
          printf("\n");
          if(outfilev1.is_open())
          {
              outfilev1<<mtx_v1[i]<<endl;    //message是程序中处理的数据
              
              
          }
      }
    ofstream outfilev2;
    outfilev2.open("/Users/h-xiao16/Desktop/数值分析与算法/mtxv2.txt"); //mtxl.txt是存放数据的文件名
    for ( unsigned i=0; i<n+3; ++i )
    {
        printf("%f",mtx_v2[i]);
        printf("\n");
        outfilev2<<mtx_v2[i]<<endl;
    }
    outfilev2.close();
    double * w1=new double [n+3];
  
    double * w2=new double [n+3];
 
 
    ofstream outfilew1;
    outfilew1.open("/Users/h-xiao16/Desktop/数值分析与算法/mtxw1.txt"); //mtxl.txt是存放数据的文件名
    ofstream outfilew2;
    outfilew2.open("/Users/h-xiao16/Desktop/数值分析与算法/mtxw2.txt"); //mtxl.txt是存放数据的文件名

   

    w1=Gause(mtx_increase,71);
    w2=Gause(mtx_increase2,71);
 //     w2=gaussin(mtx_l_change,mtx_v2_change);
//    w1=ALU(mtx_l,mtx_v1);
//    w2=ALU(mtx_l,mtx_v2);
    double **result;
    result = (double**)malloc((n+3)*sizeof(double*));
    for(int i=0;i<n+3;i++)
    {
        result[i] = (double*)malloc(2*sizeof(double));
        result[i][0]=w1[i];
        result[i][1]=w2[i];
        outfilew1<<w1[i]<<endl;
        outfilew2<<w2[i]<<endl;
    }
 outfilew1.close();
     outfilew2.close();
    return result;
}


int main()
{
    

    
    
    //加载图片路径
    Mat im_src = imread("/Users/h-xiao16/Desktop/数值分析与算法/picture/1.jpg");//特朗普
    Mat im_tg = imread("/Users/h-xiao16/Desktop/数值分析与算法/picture/2.jpg");//小孩

    Mat im_tg_change(1000,1000,CV_32FC3, Scalar(255,255,255));
   im_tg_change=im_tg;
    circle(im_tg,Point(10,100),2,CV_RGB(255,0,0),2);
    imshow("1",im_src);
    imshow("2",im_tg);
    
    FILE *fp;
    FILE *fp2;
    if((fp=fopen("/Users/h-xiao16/Desktop/数值分析与算法/picture/1.txt","rt"))==NULL)
    {
        printf("cannot open file\n");
        return 0;
    }
    for(int i=0;i<68;i++)
    {
        mypoint pt_src;
        fscanf(fp,"%f",&pt_src.x);
        fscanf(fp," ");
        fscanf(fp,"%f",&pt_src.y);
        fscanf(fp,"\n");
        src_ft.push_back(pt_src);
        
    }
    fclose(fp);
    for(int i=0;i<68;i++)
    {
        
            printf("%g",src_ft[i].x);
        printf(" ");
         printf("%g",src_ft[i].y);
        printf("\n");
    }
    if((fp2=fopen("/Users/h-xiao16/Desktop/数值分析与算法/picture/2.txt","rt"))==NULL)
    {
        printf("cannot open file\n");
        return 0;
    }
    for(int i=0;i<68;i++)
    {
        mypoint pt_src;
        fscanf(fp2,"%f",&pt_src.x);
        fscanf(fp2," ");
        fscanf(fp2,"%f",&pt_src.y);
        fscanf(fp2,"\n");
        tg_ft.push_back(pt_src);
        
    }
    fclose(fp2);
    for(int i=0;i<68;i++)
    {
        
        printf("%g",tg_ft[i].x);
        printf(" ");
        printf("%g",tg_ft[i].y);
        printf("\n");
        
        
    }

   circle(im_tg,Point(tg_ft[0].x,tg_ft[0].y),2,CV_RGB(255,0,0),2);
       circle(im_tg,Point(tg_ft[1].x,tg_ft[1].y),2,CV_RGB(255,0,0),2);
       circle(im_tg,Point(tg_ft[2].x,tg_ft[2].y),2,CV_RGB(255,0,0),2);
    imshow("3",im_tg);
    double **a=tps_result();
    double * wx=new double [n];
    double *wy=new double[n];
    double *a1=new double [2];
    double *ax=new double [2];
    double *ay= new double[2];
    for(int i=0;i<n;i++)
    {
        wx[i]=a[i][0];
        wy[i]=a[i][1];
    }
    
//    for(int i=n;i<n+3;i++)
//    {
//        cout<<a[i][0]<<endl;
//        ax[i-n]=a[i][0];
//        cout<<a[i][1]<<endl;
//        ay[i-n]=a[i][1];
//    }
//    
//
    for(int i=0;i<n;i++)
        cout<<"wx"<<wx[i]<<endl;
    for(int i=0;i<n;i++)
        cout<<"wy"<<wy[i]<<endl;
    a1[0]=a[n][0];
    a1[1]=a[n][1];
    ax[0]=a[n+1][0];
    ax[1]=a[n+1][1];
    ay[0]=a[n+2][0];
    ay[1]=a[n+2][1];
    cout<<"向量a:"<<endl;
    cout<<a1[0]<<" "<<ax[0]<<" "<<ay[0]<<endl;
        cout<<a1[1]<<" "<<ax[1]<<" "<<ay[1]<<endl;
    double dx;
    double dy;
    dx=0;dy=0;
    int tg_h=im_tg.rows;//高度600
    int tg_w=im_tg.cols;//宽度450
    vector<mypoint>tg_position;
    
    cout<<"宽度"<<tg_w<<endl<<"changdu"<<tg_h<<endl;
    for(int i=0;i<tg_w;i++)//450
    {
        for(int j=0;j<tg_h;j++)//600
        {
    im_tg_change.at<Vec3b>(j,i)[0]=0;
    im_tg_change.at<Vec3b>(j,i)[1]=0;
    im_tg_change.at<Vec3b>(j,i)[2]=0;
        }
    }
    for(int i=0;i<n;i++)
    {
        dx=0;dy=0;
    
    for(int m=0;m<n;m++)
    {
        //cout<<"距离"<<(t1-tg_ft[m]).distance()<<endl;
        dx=dx+wx[m]*tps_base_func((tg_ft[i]-tg_ft[m]).distance());
        //cout<<"wx[m]:"<<wx[m]<<endl;
        dy=dy+wy[m]*tps_base_func((tg_ft[i]-tg_ft[m]).distance());
        //cout<<"wy[m]:"<<wx[m]<<endl;
        
        
    }
        mypoint current_pt2;
        
        current_pt2.x=a1[0]+dx+tg_ft[i].x*ax[0]+tg_ft[i].y*ay[0];
        current_pt2.y=a1[1]+dy+tg_ft[i].x*ax[1]+tg_ft[i].y*ay[1];
        cout<<"特征点："<<current_pt2.x<<" "<<current_pt2.y<<endl;
    }

    for(int i=0;i<tg_w;i++)//450
    {
        for(int j=0;j<tg_h;j++)//600
        {
            mypoint t1;
            t1.x=i;
            t1.y=j;
            dx=0;
            dy=0;
        
        
            for(int m=0;m<n;m++)
            {
                //cout<<"距离"<<(t1-tg_ft[m]).distance()<<endl;
                dx=dx+wx[m]*tps_base_func((t1-tg_ft[m]).distance());
                //cout<<"wx[m]:"<<wx[m]<<endl;
                dy=dy+wy[m]*tps_base_func((t1-tg_ft[m]).distance());
                //cout<<"wy[m]:"<<wx[m]<<endl;
                

            }
           // cout<<"dx:"<<dx<<" "<<"dy:"<<dy<<endl;
           mypoint current_pt;
         
            current_pt.x=a1[0]+dx+t1.x*ax[0]+t1.y*ay[0];
            current_pt.y=a1[1]+dy+t1.x*ax[1]+t1.y*ay[1];
            tg_position.push_back(current_pt);
           // cout<<"position:"<<current_pt.x<<","<<current_pt.y<<endl;
//         if((current_pt.x<580)&&(current_pt.x>=0)&&(current_pt.y>=0)&&(current_pt.y<=800))
//         {
            int x1=floor(current_pt.x);
            int x2=ceil(current_pt.x);
            int y1=floor(current_pt.y);
            int y2=ceil(current_pt.y);
            mypoint m1;m1.x=x1;m1.y=y1;
            mypoint m2;m2.x=x1;m2.y=y2;
            mypoint n1;n1.x=x2;n1.y=y1;
            mypoint n2;n2.x=x2;n2.y=y2;
            int x_min=x1;
            int y_min=y1;
            double dis_min=min((n2-current_pt).distance(),min(min((m2-current_pt).distance(),(m1-current_pt).distance()),(n1-current_pt).distance()));
            if(dis_min==(m2-current_pt).distance())
                x_min=x1;y_min=y2;
            if(dis_min==(n1-current_pt).distance())
                x_min=x2;y_min=y1;
            if(dis_min==(n2-current_pt).distance())
                x_min=x2;y_min=y2;
//            if((m2-current_pt).distance()<(m1-current_pt).distance())
//                x_min=x1;y_min=y2;
//            else
//                if ((m3-current_pt).distance()<((m1-current_pt).distance())) {
//                    x_min=x2
//                }
         //   cout<<"像素值为"<<x_min<<" "<<y_min<<endl;
                       // cout<<"2像素值为"<<int(current_pt.x)<<" "<<int(current_pt.y)<<endl;
            Vec3b intensity=im_src.at<Vec3b>(y_min,x_min);
            
            int blue=intensity.val[0];
    
            int green=intensity.val[1];
            
            int red=intensity.val[2];
            im_tg_change.at<Vec3b>(t1.y,t1.x)[0]=blue;
            im_tg_change.at<Vec3b>(t1.y,t1.x)[1]=green;
            im_tg_change.at<Vec3b>(t1.y,t1.x)[2]=red;
//            im_tg_change.at<Vec3b>(t1.x,t1.y)[0]=blue;
//            im_tg_change.at<Vec3b>(t1.x,t1.y)[1]=green;
//            im_tg_change.at<Vec3b>(t1.x,t1.y)[2]=red;
//         }
            
//            else
//            {
//                im_tg_change.at<Vec3b>(t1.y,t1.x)[0]=0;
//                im_tg_change.at<Vec3b>(t1.y,t1.x)[1]=0;
//                im_tg_change.at<Vec3b>(t1.y,t1.x)[2]=0;
//            }
//            
        }
        

    }
    
    
    imshow("result",im_tg_change);
    cout<<"trump"<<im_src.cols<<" "<<im_src.rows;
    cout<<"over";
    waitKey();
    return 0;
    
   
    
}
