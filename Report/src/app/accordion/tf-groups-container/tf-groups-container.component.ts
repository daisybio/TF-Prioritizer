import {Component, Input, OnInit} from '@angular/core';
import {TranscriptionFactorGroup} from "../../types/types";

@Component({
  selector: 'tf-groups-container',
  templateUrl: './tf-groups-container.component.html',
  styleUrls: ['./tf-groups-container.component.css']
})
export class TfGroupsContainerComponent implements OnInit {
  @Input()
  // @ts-ignore
  tfGroups: TranscriptionFactorGroup[];

  constructor() { }

  ngOnInit(): void {
  }

}
