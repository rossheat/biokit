import 'package:meta/meta.dart';

int fibo({@required int n}) {
  Map<int, int> fiboMap = {};

  fiboMap[1] = 1;
  fiboMap[2] = 1;

  for (var i = 3; i < n + 1; i++) {
    fiboMap[i] = fiboMap[i - 1] + fiboMap[i - 2];
  }
  return fiboMap[n];
}
