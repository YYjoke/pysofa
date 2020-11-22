
[TOC]


# get-time

## #0 GitHub

```
https://github.com/Coxhuang/get_time
```

## #1 环境

```
python >= 2.7
```

## #2 使用 

### #2.1 安装

```
pip install get-time
或者
pip3 install get-time
```

### #2.1 使用


> 导入模块

```
from get_time import get_time
```

- get_time

> 使用

```
> get_time()

2019-03-31 16:44:48
```

> 参数

```
:param year: 年 (默认今年,可输入任意年份)
:param month: 月 (默认当月,可输入任意月份)
:param day: 天 (默认今天,可输入任意天数)
:param hour: 时 (默认当前,可输入任意小时)
:param minute: 分 (默认当前,可输入任意分)
:param second: 秒 (默认当前,可输入任意秒)
:param week: 星期x (默认-1,如果不等于-1,则day参数无效,可输入任意1-7)
:param last_day_of_month: 每个月的最后一天 (默认False,如果为True,则返回month最后一天)
:param type: 输出类型 (默认time,可输入"str"返回字符串类型数据)
:param detail: 是否输出时分秒? (默认输出时分秒,如果为False,则只返回 年-月-日 )
:return: time (type datetime / str)
```

- get_timestamp

> 使用


```
get_timestamp()
> 1554110680.047509
get_timestamp(detail=False)
> 1554110707
```

> 参数


```
:param detail: True输出完整的时间戳/ False输出前10位(小数点之前)
```
- get_timestamp

> 使用


```
timestamp_to_str(1554110680.047509,"%Y-%m-%d %H:%M:%S")
> 2019-04-01 17:24:40
```

> 参数


```
:param timestamp: 时间戳 
:param strformat: 转换格式 (%Y-%m-%d %H:%M:%S)
```

- str_to_timestamp

> 使用


```
str_to_timestamp("2019-04-01 17:24:40","%Y-%m-%d %H:%M:%S")
> 1554110680
```

> 参数

```
:param timestr: 时间字符串 
:param strformat: 转换格式 (%Y-%m-%d %H:%M:%S) 
```

