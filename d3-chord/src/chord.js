import {range} from "d3-array";
import {max, tau} from "./math";

function compareValue(compare) {
  return function(a, b) {
    return compare(
      a.source.value + a.target.value,
      b.source.value + b.target.value
    );
  };
}

export default function() {
  var padAngle = 0,
      sortGroups = null,
      sortSubgroups = null,
      sortChords = null;

  function chord(matrix, groupSizes) {
    var n = matrix.length,
        groupSums = [],
        groupIndex = range(n),
        subgroupIndex = [],
        chords = [],
        groups = chords.groups = new Array(n),
        subgroups = new Array(n * n),
        k,
        x,
        x0,
        dx,
        i,
        j;

    // Compute the sum.
    k = 0, i = -1; while (++i < n) {
      x = 0, j = -1; while (++j < n) {
        x += matrix[i][j];
      }
      groupSums.push(x);
      subgroupIndex.push(range(n));
      k += x;
    }

    // Sort groups…
    if (sortGroups) groupIndex.sort(function(a, b) {
      return sortGroups(groupSums[a], groupSums[b]);
    });

    // Sort subgroups…
    if (sortSubgroups) subgroupIndex.forEach(function(d, i) {
      d.sort(function(a, b) {
        return sortSubgroups(matrix[i][a], matrix[i][b]);
      });
    });

    // Convert the sum to scaling factor for [0, 2pi].
    // TODO Allow start and end angle to be specified?
    // TODO Allow padding to be specified as percentage?
    k = max(0, tau - padAngle * n) / k;
    dx = k ? padAngle : tau / n;


    // Computer the position for each group(GO term)
    x = 0, i = -1; while (++i < n) {
      x0 = x;
      x += groupSizes[di] * k
      groups[di] = {
        index: di,
        startAngle: x0,
        endAngle: x,
        value: groupSizes[di]
      };
      x += dx
    }

    // Computer the position for each sub group(GO term intersection)
    x = 0, di = -1; while (++di < n) {
      startAngle = groups[di].startAngle;
      endAngle = group[di].endAngle;
      angleRange = endAngle - startAngle;
      //calculate angle
      a = (angleRange/groupSums[di]);
      xi = startAngle, j = -1; while (++j < n){
        dj = subgroupIndex[di][j];
        v = a * matrix[di][dj];
        xj = xi + v;
        subgroups[dj * n + di] = {
          index: di,
          subindex: dj,
          startAngle: xi,
          endAngle: xj,
          value: v
        }
        xi += v;
      }

    }
    

    // Compute the start and end angle for each group and subgroup.
    // Note: Opera has a bug reordering object literal properties!
    // x = 0, i = -1; while (++i < n) {
    //   x0 = x, j = -1; while (++j < n) {
    //     var di = groupIndex[i],
    //         dj = subgroupIndex[di][j],
    //         v = matrix[di][dj],
    //         a0 = x,
    //         a1 = x += v * k;
    //     subgroups[dj * n + di] = {
    //       index: di,
    //       subindex: dj,
    //       startAngle: a0,
    //       endAngle: a1,
    //       value: v
    //     };
    //   }
    //   groups[di] = {
    //     index: di,
    //     startAngle: x0,
    //     endAngle: x,
    //     value: groupSums[di]
    //   };
    //   x += dx;
    // }

    // Generate chords for each (non-empty) subgroup-subgroup link.
    i = -1; while (++i < n) {
      j = i - 1; while (++j < n) {
        var source = subgroups[j * n + i],
            target = subgroups[i * n + j];
        if (source.value || target.value) {
          chords.push(source.value < target.value
              ? {source: target, target: source}
              : {source: source, target: target});
        }
      }
    }

    return sortChords ? chords.sort(sortChords) : chords;
  }

  chord.padAngle = function(_) {
    return arguments.length ? (padAngle = max(0, _), chord) : padAngle;
  };

  chord.sortGroups = function(_) {
    return arguments.length ? (sortGroups = _, chord) : sortGroups;
  };

  chord.sortSubgroups = function(_) {
    return arguments.length ? (sortSubgroups = _, chord) : sortSubgroups;
  };

  chord.sortChords = function(_) {
    return arguments.length ? (_ == null ? sortChords = null : (sortChords = compareValue(_))._ = _, chord) : sortChords && sortChords._;
  };

  return chord;
}
