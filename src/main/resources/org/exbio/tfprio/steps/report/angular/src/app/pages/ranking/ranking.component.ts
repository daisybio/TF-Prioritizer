import {Component, OnInit} from '@angular/core';
import {tfGroup} from "../../../interfaces";
import {DataManagerService} from "../../services/data-manager.service";
import {Subject} from "rxjs";
import {FormControl} from "@angular/forms";

@Component({
  selector: 'app-ranking',
  templateUrl: './ranking.component.html',
  styleUrls: ['./ranking.component.scss']
})
export class RankingComponent implements OnInit {
  filterControl = new FormControl('');
  tfGroups: tfGroup[] = this.dataService.groups;
  hmList: string[] = this.dataService.hms;
  hmActive: Map<string, boolean> = new Map<string, boolean>(
    this.dataService.hms.map(hm => [hm, true])
  );
  tfGroupRanking: Subject<tfGroup[]> = new Subject<tfGroup[]>();
  openGroup: Subject<tfGroup> = new Subject<tfGroup>();
  filteredTfGroups: Subject<Map<tfGroup, boolean>> = new Subject<Map<tfGroup, boolean>>();

  constructor(private dataService: DataManagerService) {
  }

  ngOnInit(): void {
    this.filterControl.valueChanges.subscribe(value => {
      this.updateVisibleGroups(value);
    });
  }

  ngAfterViewInit(): void {
    setTimeout(() => {
      this.recalculateScores();
    }, 0);

    setTimeout(() => {
      this.updateVisibleGroups(null);
    }, 0);
  }


  async recalculateScores() {
    let groupDcg: Map<tfGroup, number> = this.tfGroups.map(group => {
      let dcgSum = this.hmList.filter(hm => this.hmActive.get(hm)).map(hm => group.hmDcg[hm])
        .filter(dcg => dcg >= 0).reduce((a, b) => a + b, 0);
      return {"group": group, "dcgSum": dcgSum}
    }).reduce((map, obj) => {
      map.set(obj.group, obj.dcgSum);
      return map;
    }, new Map<tfGroup, number>());

    let groupDcgSorted = Array.from(groupDcg.entries()).sort((a, b) => b[1] - a[1]);
    let groupDcgSortedNames = groupDcgSorted.map(entry => entry[0]);

    this.tfGroupRanking.next(groupDcgSortedNames);
  }

  private async updateVisibleGroups(filter: string | null) {
    this.filteredTfGroups.next(this.tfGroups.reduce((map, group) => {
      map.set(group, !filter || group.symbol.toLowerCase().includes(filter.toLowerCase()));
      return map;
    }, new Map<tfGroup, boolean>()));
  }
}
