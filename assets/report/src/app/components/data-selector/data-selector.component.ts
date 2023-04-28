import {Component, Input} from '@angular/core';
import {Subject} from "rxjs";

@Component({
  selector: 'app-data-selector',
  templateUrl: './data-selector.component.html',
  styleUrls: ['./data-selector.component.scss']
})
export class DataSelectorComponent {
  @Input()
  data: {} = {};
  @Input()
  activeData!: Subject<any>;
  @Input()
  maxDepth?: number;
  levelExisting: string[][] = [];
  actives: Subject<string[]> = new Subject<string[]>();
  currentActives: string[] = [];
  levelAllowed: Subject<Map<string, boolean>[]> = new Subject<Map<string, boolean>[]>();
  imagePath: Subject<string> = new Subject<string>();
  dropdown: boolean = false;

  ngOnInit() {
    this.levelExisting = this.extractKeys(this.data);
    if (this.maxDepth) {
      this.levelExisting = this.levelExisting.slice(0, this.maxDepth);
    }

    this.dropdown = this.levelExisting[this.levelExisting.length - 1].length > 10;

    this.actives.subscribe(actives => {
      let allowed: Map<string, boolean>[] = [];
      for (let i = 0; i < actives.length; i++) {
        allowed.push(new Map<string, boolean>());
        let found: string[] = this.getChildKeys(this.data, actives.slice(0, i));
        this.levelExisting[i].forEach(key => allowed[i].set(key, found.includes(key)));
      }

      this.levelAllowed.next(allowed);
    });

    this.actives.subscribe(actives => {
      this.activeData.next(this.getData(this.data, actives));
    })
  }

  ngAfterViewInit(): void {
    setTimeout(() => {
      let actives: string[] = [];
      for (let i = 0; i < this.levelExisting.length; i++) {
        actives.push(this.getChildKeys(this.data, actives).sort()[0]);
      }

      this.currentActives = actives;

      this.actives.next(actives);
    }, 0);
  }

  setLevelActive(level: number, active: string) {
    if (active === null || active === undefined || active === "") {
      this.actives.next(this.currentActives);
      return;
    }

    console.log("setLevelActive", level, active);
    this.currentActives[level] = active;

    for (let i = level + 1; i < this.levelExisting.length; i++) {
      let levelActive = this.currentActives[i];
      let childKeys = this.getChildKeys(this.data, this.currentActives.slice(0, i));
      if (childKeys.includes(levelActive)) {
        continue;
      }
      this.currentActives[i] = childKeys.sort()[0];
    }

    this.actives.next(this.currentActives);
  }


  private extractKeys(obj: any, level: number = 0, result: string[][] = []): string[][] {
    if (!result[level]) result[level] = []
    for (const key in obj) {
      if (obj.hasOwnProperty(key)) {
        if (!result[level].includes(key))
          result[level].push(key)
        if (typeof obj[key] === "object") {
          this.extractKeys(obj[key], level + 1, result);
        }
      }
    }
    result.forEach(level => level.sort())
    return result;
  }

  private getChildKeys(obj: any, keys: string[]): string[] {
    let currentValue = obj;
    for (const key of keys) {
      currentValue = currentValue[key];
      if (currentValue === undefined) {
        return [];
      }
    }
    return Object.keys(currentValue);
  }

  private getData(obj: any, keys: string[]): any {
    let currentValue = obj;
    for (const key of keys) {
      currentValue = currentValue[key];
      if (currentValue === undefined) {
        return "";
      }
    }
    return currentValue;
  }
}
