import {Component, Input, OnInit} from '@angular/core';

@Component({
  selector: 'data-entry',
  templateUrl: './data-entry.component.html',
  styleUrls: ['./data-entry.component.css']
})
export class DataEntryComponent implements OnInit {
  @Input()
  // @ts-ignore
  title: string;

  @Input()
  // @ts-ignore
  information: string;

  informationVisible: boolean = false;

  constructor() { }

  ngOnInit(): void {
  }

  toggleInformationVisibility() {
    this.informationVisible = !this.informationVisible;
  }
}
