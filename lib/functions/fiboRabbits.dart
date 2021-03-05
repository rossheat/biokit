import 'package:meta/meta.dart';

fiboRabbits({@required int months, @required int pairs}) {
  Map<int, int> fiboMap = {};

  fiboMap[1] = 1;
  fiboMap[2] = 1;

  for (var month = 3; month <= months; month++) {
    fiboMap[month] = (fiboMap[month - 2] * pairs) + fiboMap[month - 1];
  }

  return fiboMap[months];
}
