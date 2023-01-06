#' The Use of Marginal Distributions in Conditional Forecasting
#' @description A new way to predict time series using the marginal distribution table in the absence of the significance of traditional models.
#' @return the output from \code{ff()}
#' @export
#' @param dt data frame
#' @param m the number of time series
#' @param w the number of predicted values
#' @param q1 matrix independent time series values #In the case of m=2, enter the independent string values as follows(matrix(c())),In the case of m=3, enter the independent string values as follows(matrix(c(),w,m-1,byrow=T))
#' @param n number of values
#' @usage ff(dt,m,w,n,q1)
#' @examples
#' x=rnorm(17,10,1)
#' y=rnorm(17,10,1)
#' data=data.frame(x,y)
#' print("Enter independent time series values")
#' q1=list(q=matrix(c(scan(,,quiet=TRUE)),1,2-1))
#' 10.5
#'
#'
#' ff(data,2,1,17,q1)
ff=function(dt,m,w,n,q1)
    { qq=list(Number_of_categories= ceiling (2.5*(n^(0.25))),mx=matrix(nrow=m,ncol=1),mn= matrix(nrow=m,ncol=1),rn= matrix(nrow=m,ncol=1),l= matrix(nrow=m,ncol=1)
       ,xd=matrix(nrow=ceiling (2.5*(n^(0.25))),ncol=m),xu= matrix(nrow=ceiling (2.5*(n^(0.25))),ncol=m),xx=matrix(nrow=ceiling (2.5*(n^(0.25))),ncol=m),r=matrix(nrow=ceiling (2.5*(n^(0.25))),ncol=ceiling (2.5*(n^(0.25)))^(m-1)),
         rr=matrix(nrow=1,ncol=ceiling (2.5*(n^(0.25)))^(m-1)),rr3=matrix(nrow=1,ncol=ceiling (2.5*(n^(0.25)))^(m-1)),predicted_values=matrix(nrow = w,ncol = 1))

  n=nrow(dt)
  if (is.null(n) || is.na(n) || is.infinite(n)) {
    stop(terminatepg)
  }
  terminatepg = "Terminating program ... \n"
  origmsg = "Here is the original message: "

  dt = tryCatch(
    {
      as.data.frame(dt)
    },
    warning = function(warning_msg) {
      message("dt The entered value is non-numeric, please try again")
      message(origmsg)
      message(paste(warning_msg, "\n"))
      return(NULL)
    }
  )
  m = tryCatch(
    {
      as.integer(m)
    },
    warning = function(warning_msg) {
      message("m The entered value is non-numeric, please try again")
      message(origmsg)
      message(paste(warning_msg, "\n"))
      return(NULL)
    }
  )
w = tryCatch(
    {
      as.integer(w)
    },
    warning = function(warning_msg) {
      message("w The entered value is non-numeric, please try again")
      message(origmsg)
      message(paste(warning_msg, "\n"))
      return(NULL)
    }
    )

  nameMsg = "Number of categories format not supported! "

       for(i in 1:m)
    {
      qq$mx[i]=max(dt[,i])
      qq$mn[i]=min(dt[,i])
      qq$rn[i]=qq$mx[i]-qq$mn[i]
    }
    if (is.null(qq$mx[i]) || is.na(qq$mx[i]) || is.infinite(qq$mx[i])) {
      stop(terminatepg)
    }
    if (is.null(qq$mn[i]) || is.na(qq$mn[i]) || is.infinite(qq$mn[i])) {
      stop(terminatepg)
    }
    if (is.null(qq$rn[i]) || is.na(qq$rn[i]) || is.infinite(qq$rn[i])) {
      stop(terminatepg)
    }

    for(i in 1:m)
    {
    qq$l[i]=(qq$rn[i]/(qq$Number_of_categories))
    }
    if (is.null(qq$l[i]) || is.na(qq$l[i]) || is.infinite(qq$l[i])) {
      stop(terminatepg)
    }

    for(i in 1:m)
    {
      qq$xd[1,i]=qq$mn[i]
      qq$xu[1,i]=qq$mn[i]+qq$l[i]
      qq$xx[1,i]=(qq$xd[1,i]+qq$xu[1,i])/2
      for(j in 2:qq$Number_of_categories)
      {
        qq$xd[j,i]=qq$xd[j-1,i]+qq$l[i]
        qq$xu[j,i]=qq$xu[j-1,i]+qq$l[i]
        qq$xx[j,i]=(qq$xd[j,i]+qq$xu[j,i])/2
      }}
      for(i in 1:m)
    {qq$xu[qq$Number_of_categories,i]=qq$xu[qq$Number_of_categories,i]+0.5}
    if (is.null(qq$xd[i]) || is.na(qq$xd[i]) || is.infinite(qq$xd[i])) {
      stop(terminatepg)
    }
    if (is.null(qq$xu[i]) || is.na(qq$xu[i]) || is.infinite(qq$xu[i])) {
      stop(terminatepg)
    }
    if (is.null(qq$xx[i]) || is.na(qq$xx[i]) || is.infinite(qq$xx[i])) {
      stop(terminatepg)
    }

    if(m==2)
    {
      for(i in 1:qq$Number_of_categories)
      {
        for(s in 1:qq$Number_of_categories)
        {
          qq$r[i,s]=0
          for(j in 1:n)
          {
            ss1= dt[j,1]>=qq$xd[s,1] & dt[j,1]<qq$xu[s,1]
            ss2= dt[j,2]>=qq$xd[i,2] & dt[j,2]<qq$xu[i,2]
            ss=(ss1 & ss2)
            if(ss)
            {qq$r[i,s]=qq$r[i,s]+1}
          }}}
      qq$r=round(qq$r/n,2)
      qq$rr=apply(qq$r,2,sum)
      qq$rr3=apply(qq$r,1,sum)
      if (is.null(qq$r[i]) || is.na(qq$r[i]) || is.infinite(qq$r[i])) {
        stop(terminatepg)
      }
      if (is.null(qq$rr[i]) || is.na(qq$rr[i]) || is.infinite(qq$rr[i])) {
        stop(terminatepg)
      }
      if (is.null(qq$rr3[i]) || is.na(qq$rr3[i]) || is.infinite(qq$rr3[i])) {
        stop(terminatepg)
      }
      for(i in 1:w)
      {
      y=0
      sss=0
      for(s in 1:qq$Number_of_categories)
      {
        sss=sss+1
        if(q1[i]>=qq$xd[s,1] & q1[i]<qq$xu[s,1])
        {if(qq$rr[sss]!=0)
        {
          y=y+(qq$r[,sss]*(qq$xx[,2]))/qq$rr[sss]}
          else
          {for(sss in 1:qq$Number_of_categories)
            y=y+(qq$rr3[sss]*(qq$xx[sss,2]))}}
        qq$predicted_values[i]=sum(y)
      }}
      return(qq)
      if (is.null(qq$predicted_values[i]) || is.na(qq$predicted_values[i]) || is.infinite(qq$predicted_values[i])) {
        stop(terminatepg)

}}
  else
    {

      sss=0
      for(i in 1:qq$Number_of_categories)
      {
        for(s in 1:qq$Number_of_categories)
        {
          sss=sss+1
          for(h in 1:qq$Number_of_categories)
          {
            qq$r[h,sss]=0
            for(j in 1:n)
            {
              ss1= (dt[j,1]>=qq$xd[i,1] & dt[j,1]<qq$xu[i,1])
              ss2= (dt[j,2]>=qq$xd[s,2] & dt[j,2]<qq$xu[s,2])
              ss3= (dt[j,3]>=qq$xd[h,3] & dt[j,3]<qq$xu[h,3])
              ss=(ss1 & ss2 & ss3)
              if(ss)
              {qq$r[h,sss]=qq$r[h,sss]+1}
            }}}}
      qq$r=round(qq$r/n,2)
      qq$rr=apply(qq$r,2,sum)
      qq$rr3=apply(qq$r,1,sum)

      for(i in 1:w)
      {
        y=0
        sss=0
        for(s in 1:qq$Number_of_categories)
        {
          for(h in 1:qq$Number_of_categories)
          {
            sss=sss+1
            if(q1[i,1]>=qq$xd[s,1] & q1[i,1]<qq$xu[s,1])
            {
              if(q1[i,2]>=qq$xd[h,2] & q1[i,2]<qq$xu[h,2])
              {
                if(qq$rr[sss]!=0)
                {
                  y=y+(qq$r[,sss]*qq$xx[,3])/qq$rr[sss]}
                else
                {for(sss in 1:qq$Number_of_categories)
                  y=y+(qq$rr3[sss]*qq$xx[sss,3])}}
              qq$predicted_values[i]=sum(y)
            }}}}
        return(qq)
        if (is.null(qq$predicted_values[i]) || is.na(qq$predicted_values[i]) || is.infinite(qq$predicted_values[i])) {
          stop(terminatepg)}
        }}
